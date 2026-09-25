# BOSS DSL for 2401.08252v3: First observation of psi(3686) → Omega- K+ anti-Xi0 + c.c.
# Data: (27.12±0.14)×10^8 psi(3686) events at 3.686 GeV
# Uses partial reconstruction: reconstruct Omega- K+, infer anti-Xi0 via recoil mass
# Omega- → Lambda K-, Lambda → p pi-

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(3686) → Omega- K+ anti-Xi0
# Omega- → Lambda K-, Lambda → p pi-
# anti-Xi0 is NOT reconstructed (partial_rec)
# Need recID mapping for partial_rec:
# RecID list:
#   0 => psi(3686)    (top mother — skip)
#   1 => Omega-
#   2 => K+
#   3 => anti-Xi0      (untagged recoil — skip)
#   4 => Lambda
#   5 => K-            (from Omega-)
#   6 => p             (from Lambda)
#   7 => pi-           (from Lambda)
# Then reconstruct: Omega-, K+, Lambda, K-, p, pi-
# Miss: anti-Xi0 (recID 3)
# recIDs to reconstruct: [1, 2, 4, 5, 6, 7]  (everything except 0 and 3)
# Or use partial_miss([3]) — easier

decay_card = <<~DECAYCARD
  Decay psi(3686)
  1.0000 Omega- K+ anti-Xi0 PHSP;
  Enddecay
  Decay Omega-
  1.0000 Lambda K- VSS;
  Enddecay
  Decay Lambda
  1.0000 p+ pi- VSS;
  Enddecay
  Decay anti-Xi0
  1.0000 anti-Lambda pi0 PHSP;
  Enddecay
  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_omegak_xibar0"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# Event selection: reconstruct Omega- → Lambda(p pi-) K-
# anti-Xi0 by recoil mass RM(Omega- K+)
event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"   # p, K+ at minimum
    nChrn     ">=2"   # pi-, K- at minimum (from Omega-)
  end
  .pid(method: :probability) do
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
  end
  .pid(method: :probability) do
    identify :kaon, against: [:pion]
  end
  .pid(method: :probability) do
    identify :pion, against: [:kaon, :proton]
  end
  # Lambda vertex fit: p pi- pair
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Now do partial reconstruction:
  # recIDs: 0=psi(3686), 1=Omega-, 2=K+, 3=anti-Xi0, 4=Lambda, 5=K-, 6=p, 7=pi-
  # Reconstruct: Omega-(1), K+(2), Lambda(4), K-(5), p(6), pi-(7)
  # Miss: anti-Xi0(3)
  .partial_miss([3]) do
    require_recoil_mass 1.282, 1.352
    # Also require best combination by Omega- mass
    best_combination_by_mass :Omega, 1.67245
  end

algorithm = Algorithm.new("PsipToOmegaKXibar0")
algorithm
  .set_header(["PsipToOmegaKXibar0/PsipToOmegaKXibar0.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .note(:kplus_track_quality,
    "K+ from IP must have |Vz|<10cm and |Vxy|<1cm. If multiple K+, select highest CL_K. " \
    "Not expressible in DSL — need custom track selection.")
  .note(:lambda_mass_window,
    "Lambda mass window: M(p pi-) in [1.111, 1.121] GeV/c^2. Applied in ROOT stage.")
  .note(:omega_mass_window,
    "Omega- reduced mass M_Omega = M(Lambda K-) - M(p pi-) + M_Lambda_PDG " \
    "in [1.663, 1.681] GeV/c^2 (~6 sigma). Signal/sideband regions applied in ROOT.")
  .note(:omega_vertex_fit,
    "Omega- vertex fit chi2 < 200. Omega- reconstructed from Lambda + K- with common vertex.")
  .note(:best_omega_selection,
    "If multiple Omega- candidates, choose minimum |M_Omega - M_Omega_PDG|. " \
    "Multiple candidate rate ~1.7%, with <4% wrong selection.")
  .note(:kplus_pid_condition,
    "K+ PID: CL_K > CL_pi and CL_K > 0. Proton PID: CL_p > CL_K, CL_p > CL_pi, CL_p > 0.001. " \
    "Remaining tracks assigned as pions.")
  .note(:continuum_subtraction,
    "Continuum background from 3.773 GeV data: scale factor f_c applied with " \
    "luminosity and cross-section normalization. n=1 assumed with 1/s dependence.")
  .note(:no_intermediate_states,
    "No evidence for intermediate resonances (Xi(2250)^0 → Omega- K+). " \
    "PHSP signal MC shows acceptable agreement with data.")
  .with_decay_card(decay_card)
  .apply(event_selection)
  .execute_on([psip_data, psip_incMC, exMC])