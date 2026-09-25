# Measurement of χ_cJ → φ φ η branching fractions
# arXiv:1911.02988v1 — ψ(3686) radiative decays
# (448.1 ± 2.9) × 10⁶ ψ(3686) events

### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

### Decay card for ψ(3686) → γ χ_cJ, χ_cJ → φ φ η, φ → K+K-, η → γγ ###
decay_card = <<~DECAYCARD
  Decay psi(3686)
  1.0000 gamma chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.0000 phi phi eta PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC for signal ###
exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_psip_gamma_chic0_phiphi_eta"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("ChicJPhiPhiEta")
alg.set_header(["ChicJPhiPhiEtaAlg/ChicJPhiPhiEta.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"      # 2 K+ from two φ decays
  nChrn "==2"      # 2 K- from two φ decays
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=3"       # 1 radiative photon + 2 from η → γγ
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==2"
  nkm "==2"
end
# Assign charged tracks — all are identified as kaons
.remove([:prp <= :chrgp, :prm <= :chrgn])
.assign({chrgp: :kp, chrgn: :km})
# Reconstruct η → γγ via Kalman fit
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end
# 4C kinematic fit: γ_radiative + K+ + K- + K+ + K- + η
# Best photon assigned as radiative, remaining γγ → η already reconstructed
.kinematic_fit([:gamma, :kp, :km, :kp, :km, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

# Attach notes for BOSS-side procedures not expressible in formal DSL constructs
alg.note(:eta_signal_window, "η signal region: M(γγ) ∈ [0.52, 0.58] GeV/c²; best γγ pair chosen as closest to nominal η mass")
   .note(:phi_signal_window, "φ signal region: M(K+K-) ∈ [1.005, 1.035] GeV/c²; φ candidates chosen by minimising ΔM² = (M1-m_φ)² + (M2-m_φ)²")
   .note(:radiative_photon_selection, "the photon not used in η reconstruction is treated as the radiative photon from ψ(3686) → γ χ_cJ")
   .note(:jpsi_eta_veto, "recoil mass against η required < 3.05 GeV/c² to suppress ψ(3686) → η J/ψ, J/ψ → γ φ φ")
   .note(:pi0_veto, "all γγ combinations vetoed in [0.115, 0.150] GeV/c² to suppress π⁰ background")
   .note(:phi_eta_veto, "M(γ η) vetoed in [1.00, 1.04] GeV/c² to suppress ψ(3686) → φ φ φ where φ → γ η")
   .note(:helix_correction, "helix parameter correction applied to charged tracks from KS0 control sample")
   .note(:e1_weighting, "E1 transition weighting factor (E_γ/E_γ0)³ applied to M(φ K+K- γγ) spectra for detection efficiency correction")
   .note(:chic_signal_regions, "χ_c0: M(φφη) ∈ [3.38,3.45]; χ_c1: [3.48,3.54]; χ_c2: [3.54,3.60] GeV/c²")

alg.with_decay_card(decay_card).apply(event_selection)
root_files = alg.execute_on([psip_data, psip_incMC, exMC])