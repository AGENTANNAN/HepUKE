# =============================================================================
# e+e- -> gamma omega J/psi, omega -> pi+ pi- pi0, pi0 -> gamma gamma,
# J/psi -> l+ l- (l = e, mu);  energy scan 4.009 - 4.600 GeV (19 points)
# =============================================================================

### ---------------------------- Dataset preparation --------------------------- ###
# 19 scan points between 4.009 and 4.600 GeV (BOSS release 703)
scan_energies = [4009, 4180, 4190, 4200, 4210, 4220, 4230, 4246, 4260,
                 4270, 4280, 4310, 4360, 4390, 4420, 4470, 4530, 4575, 4600]
scan_data  = scan_energies.map { |e| DatasetManager.real_data.find("703_#{e}") }    # real data
scan_incMC = scan_energies.map { |e| DatasetManager.inclusive_mc.find("703_#{e}") } # inclusive MC

# Signal decay card (KKMC; top mother = virtual photon)
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 gamma omega J/psi      PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0            OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma            PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e-                  PHOTOS VLL;
    0.500 mu+ mu-                PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# One exclusive signal MC sample per scan point (same card, cross section and statistics)
exMCs = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
    config.sample_name   = "exmc_scan_gamma_omega_jpsi"
    config.events        = 100000
    config.decay_card    = decay_card_signal
    config.cross_section = :default
end

### ----------------------------- Event selection (BOSS) ----------------------- ###
alg_name = "GammaOmegaJpsi"
gamma_omega_jpsi = Algorithm.new(alg_name)
gamma_omega_jpsi.set_header(["#{alg_name}Alg/#{alg_name}.h"])
                .set_constant({"ECMS" => [:double, 4.26]})
                .set_alias({"std::vector<double>" => "Vdouble"})
                .note(:beam_energy_scan,
                      "the analysis runs over 19 scan points from 4.009 to 4.600 GeV; the beam " \
                      "energy (and therefore the CMS four-momentum entering the kinematic fit) " \
                      "differs at every point and is read from the per-run conditions/database " \
                      "rather than the single ECMS constant used in this algorithm")
                .note(:pid_correction_method,
                      "momentum-based PID: tracks with p > 1 GeV/c are treated as J/psi lepton " \
                      "candidates and tracks with p < 1 GeV/c as omega pion candidates; e/mu " \
                      "separation uses EMC energy (muon E_EMC < 0.35 GeV, electron E_EMC > 1.1 GeV). " \
                      "identify_high_momentum_leptons only exposes a single EMC threshold " \
                      "(set to 1.1 GeV), so the muon E_EMC < 0.35 GeV requirement is enforced " \
                      "explicitly in the BOSS algorithm")

event_selection = Selection.new
event_selection
    # ---- charged track selection ----
    .select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93
        Vr        1.0       # Vr < 1 cm
        Vz        10.0      # |Vz| < 10 cm
        nChrp     "==2"     # pi+ and l+
        nChrn     "==2"     # pi- and l-
        nNet      "==0"     # net charge zero
    }
    # ---- photon selection ----
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025     # 25 MeV, barrel
        energyThreshold_e 0.050     # 50 MeV, end-cap
        angle_to_track    20.0      # > 20 deg from any charged track
        nGam              ">=3"     # at least 3 photons
        nGam              "<=20"    # at most 20 photons
    }
    # ---- momentum-based lepton / pion separation ----
    .pid(method: :probability) {
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,   # p > 1 GeV/c -> l
                                        treat_as_electron_if_energy_above: 1.1   # E_EMC > 1.1 GeV -> e
        identify :pion, against: [:kaon, :proton]                                # p < 1 GeV/c -> pi
        nlp  "==1"      # one l+
        nlm  "==1"      # one l-
        npip "==1"      # one pi+
        npim "==1"      # one pi-
    }
    # ---- pi0 reconstruction: two photons mass-constrained to m(pi0) ----
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0     ">=1"
    }
    # ---- 5C kinematic fit: four-momentum conservation + pi0 mass constraint ----
    # the pi0 mass constraint is carried into the fit by the pre-reconstructed pi0
    .kinematic_fit([:gamma, :pi0, :pip, :pim, :lp, :lm]) {
        nominal                 # best combination (smallest chi2) fixes the radiative photon
        constrain_four_momentum
        chi2_cut 100
    }

gamma_omega_jpsi.with_decay_card(decay_card_signal).apply(event_selection)
root_files = gamma_omega_jpsi.execute_on(scan_data + scan_incMC + exMCs)