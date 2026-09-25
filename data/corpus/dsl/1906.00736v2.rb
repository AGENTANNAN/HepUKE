# BESIII DSL: Study of ψ(3686) → φ π+π− η with φ → K+K−, η → γγ
# ArXiv: 1906.00736v2
# Also measures φη', φf1(1285), φη(1405)
# ψ(3686) at 3.686 GeV

### Dataset preparation ###

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for ψ(3686) → φ π+π− η, φ → K+K−, η → γγ
decay_card_phi_pipieta = <<~DECAYCARD
  Decay psi(2S)
  1.000 phi pi+ pi- eta PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for signal
sig_phi_pipieta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_phi_pipieta_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_phi_pipieta
  config.cross_section   = :default
end

### Event selection (BOSS) ###

alg_phi_pipieta = Algorithm.new("PhiPiPiEta")
alg_phi_pipieta.set_header(["PhiPiPiEtaAlg/PhiPiPiEta.h"])
                .set_constant({ "ECMS" => [:double, 3.686] })

sel_phi_pipieta = Selection.new
  # Select charged tracks: 2 positively charged + 2 negatively charged
  .select_track {
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
  }
  # Select photons: at least 2 photons
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  }
  # PID: identify kaons and pions using probability method
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kp, against: [:pion, :proton]
    identify :km, against: [:pion, :proton]
    identify :pip, against: [:kaon, :proton]
    identify :pim, against: [:kaon, :proton]
    nkp "==1"
    nkm "==1"
    npip "==1"
    npim "==1"
  }
  # Remove identified particles from generic lists
  .remove([:kp <= :chrgp, :km <= :chrgn, :pip <= :chrgp, :pim <= :chrgn])
  # 4C kinematic fit: ψ(3686) → K+K−π+π−γγ
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 80
  }

alg_phi_pipieta
  .note(:phi_mass_window, "phi signal region: |M(K+K-) - m_phi| < 0.015 GeV/c^2")
  .note(:eta_mass_window, "eta signal region: |M(γγ) - m_eta| < 0.040 GeV/c^2")
  .note(:jpsi_veto, "J/psi veto: |M(K+K-γγ) - M_Jpsi| > 0.05 GeV/c^2 and |M(K+K-π+π-) - M_Jpsi| > 0.05 GeV/c^2 when M(π+π-η) > 1.1 GeV/c^2; suppresses ψ(3686)→π+π-J/ψ background")
  .note(:continuum_subtraction, "continuum e+e-→φπ+π-η background estimated from ψ(3770) data at 3.773 GeV and subtracted using scale factor based on luminosity, cross section and efficiency ratios")
  .note(:sideband_background, "η-φ 2D sidebands: η sidebands [0.448,0.488] or [0.608,0.648] GeV/c^2, φ sideband [1.045,1.075] GeV/c^2; used for background estimation")
  .note(:helicity_amplitude_model, "ψ(3686)→φη' and ψ(3686)→φη(1405) generated with helicity amplitude model; ψ(3686)→φf1(1285) uses same model as J/ψ→φf1(1285)")
  .note(:unbinned_ml_fit, "Unbinned maximum likelihood fit over M(π+π-η): range [0.85,1.10] GeV/c^2 for η', and [1.1,2.2] GeV/c^2 for f1(1285) + η(1405) with Breit-Wigner convolved with Gaussian resolution")
  .with_decay_card(decay_card_phi_pipieta)
  .apply(sel_phi_pipieta)

alg_phi_pipieta.execute_on([psip_data, psip_incMC, sig_phi_pipieta])