# Paper 2308.15206v3: ψ(3686)→K-ΛΞ++c.c. partial wave analysis of Ξ* hyperons
# Partial reconstruction: Λ from recoil mass RM(K-Ξ+)
# Single energy: ψ(3686), 448.1M events, BOSS 709

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay chain: ψ(3686) → K- Λ anti-Ξ+ (c.c.)
#   anti-Ξ+ → anti-Λ π+
#     anti-Λ → anti-p π+
#   Λ → p π-  (missing — inferred from recoil)

decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.0000 K- Lambda anti-Xi+  PHSP;
  Enddecay
  Decay anti-Xi+
  1.0000 anti-Lambda pi+  PHSP;
  Enddecay
  Decay anti-Lambda
  1.0000 anti-p- pi+  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  End
DECAYCARD

algorithm = Algorithm.new("Psip_K_Lambda_Xi", "00-00-01")
  .set_header(["Psip_K_Lambda_Xi/Psip_K_Lambda_Xi.h"])
  .set_constant(ECMS: 3.686109)
  .with_decay_card(decay_card)

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_K_Lambda_Xi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      ">=4"
    nChrp     ">=2"
    nChrn     ">=2"
  end
  .pid(method: :probability) do
    identify :proton, against: [:kaon, :pion]
    identify :kaon, against: [:pion]
    nprm "==1"
    nkm  "==1"
  end
  .remove([:prm <= :chrgn, :km <= :chrgn])
  .pid(method: :probability) do
    identify :pion, against: [:kaon]
    npip "==2"
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .build_virtual_particle(:Xi_plus_bar, from: [:Lambda_bar, :pip])
  .kinematic_fit([:km, :Xi_plus_bar]) do
    miss_track_of :Lambda
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algorithm
  .note(:loose_track_cuts_xidaughters, "Xi+ daughter tracks use looser Vz (15 cm) and Vr (10 cm) cuts in addition to standard MDC track selection; applied via BOSS-side track-parameter override")
  .note(:lambda_mass_window, "anti-Lambda candidates selected with M(pbar pi+) in [1.110, 1.121] GeV/c^2")
  .note(:xi_decay_length, "anti-Xi+ decay length > 0.5 cm required")
  .note(:xi_mass_window, "M(Lambda_bar pi+) in [1.315, 1.330] GeV/c^2 for anti-Xi+ candidates")
  .note(:recoil_mass_window, "RM(K- anti-Xi+) in [1.080, 1.140] GeV/c^2 to select prompt Lambda candidates")
  .apply(selection)

algorithm.execute_on([psip_data, psip_incMC, exMC])