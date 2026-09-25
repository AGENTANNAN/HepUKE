### dataset description ###
psi3686_data = DatasetManager.real_data.find("709_3686")
psi3686_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay: psi(3686) -> gamma chi_c0, chi_c0 -> omega phi, omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
decay_card = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma chi_c0    VSP_PWAVE;
    Enddecay

    Decay chi_c0
    1.0000 omega phi    HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0    PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_chicJ_omega_phi"
    config.related_dataset = psi3686_data
    config.events = 100000
    config.decay_card = decay_card
    config.cross_section = :default
end

### chi_cJ (J=0,1,2) -> omega phi ###
alg = Algorithm.new("ChicJOmegaPhi")
alg.set_header(["ChicJOmegaPhiAlg/ChicJOmegaPhi.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new
sel.select_track {
        cos_theta  0.93
        Vz         10.0
        Vr         1.0
        nChrp      "==2"
        nChrn      "==2"
        nNet       "==0"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=3"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        nkp ">=0"; nkm ">=0"
    }
    .remove([:kp <= :chrgp, :km <= :chrgn])
    .assign({chrgp: :pip, chrgn: :pim})
    .kalman_kinematic_fit([:gamma, :gamma]) {
        build_virtual_particle(:pi0).by_minimizing_mass_difference
    }
    .kinematic_fit([:gamma, :pi0, :kp, :km, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 60
    }

alg
  .note(:kaon_pid, "at least one charged track identified as kaon via dE/dx+TOF; second kaon found by minimizing |M(K+K-)-M_phi| among all identified kaons and opposite-charge tracks; remaining two tracks assigned as pions; DSL approximates by identifying kaons with probability method and removing them from charged track lists")
  .note(:pi0_selection, "pi0 candidate selected from three gamma-gamma combinations choosing pair with minimum |M(gamma gamma) - M_pi0|; applied to 4C-fit-corrected momenta; DSL approximates via kalman_kinematic_fit before 4C fit")
  .note(:competing_hypothesis_veto, "chi2_4C(3gamma K+K-pi+pi-) < chi2_4C(4gamma K+K-pi+pi-) applied to suppress background with extra photon in final state; 4gamma hypothesis tested when >3 photon candidates present")
  .note(:jpsi_veto, "pi+pi- recoil mass veto: |M_recoil(pi+pi-) - M_J/psi| > 8 MeV/c^2 to suppress psi(3686) -> pi+pi- J/psi background")
  .note(:omega_mass_window, "omega candidate with |M(pi+pi-pi0) - M_omega| < 0.05 GeV/c^2 after 4C kinematic fit; applied to 4C-fit-corrected momenta")
  .note(:phi_mass_window, "phi candidate with |M(K+K-) - M_phi| < 0.015 GeV/c^2 after 4C kinematic fit; applied to 4C-fit-corrected momenta")
  .note(:sideband_method, "signal region D defined by omega and phi mass windows; five sideband regions A, B, C defined in 2D M(K+K-) vs M(pi+pi-pi0) plane; peaking backgrounds (chi_cJ->omega K+K-, phi pi+pi-pi0, non-resonant K+K-pi+pi-pi0) estimated from sidebands with MC-derived scaling factors")
  .note(:signal_extraction, "chi_cJ yields extracted via simultaneous unbinned ML fit to M(K+K-pi+pi-pi0) distributions in signal and sideband regions; signal shape = MC simulated shape convoluted with Gaussian for data-MC resolution difference; background = polynomial + sideband-normalized peaking components")
  .note(:helamp_model, "signal MC: psi(3686)->gamma chi_cJ generated as E1 transition (VSP_PWAVE); chi_cJ->omega phi generated with HELAMP model using same helicity amplitudes as chi_cJ->phi phi; omega and phi decays via PHSP")
  .note(:chi_cJ_states, "three chi_cJ states (J=0,1,2) share identical final state and selection; analyzed simultaneously via fit to M(K+K-pi+pi-pi0); chi_c1 observed at 12.3 sigma, chi_c2 at 4.8 sigma, chi_c0 measured with improved precision")
  .note(:continuum, "continuum background evaluated using sqrt(s)=3.65 GeV data (~1/15 of psi(3686) luminosity); found to be negligible after selection")
  .with_decay_card(decay_card)
  .apply(sel)

alg.execute_on([psi3686_data, psi3686_incMC, exMC])