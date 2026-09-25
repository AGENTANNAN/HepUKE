# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
# 19 real-data energy points covering sqrt(s) = 4.310 - 4.946 GeV (BOSS 703/706/707)
data_points = [
  DatasetManager.real_data.find("703_4360"),  # ~4.358 GeV
  DatasetManager.real_data.find("703_4390"),  # ~4.387 GeV
  DatasetManager.real_data.find("703_4420"),  # ~4.416 GeV
  DatasetManager.real_data.find("703_4470"),  # ~4.467 GeV
  DatasetManager.real_data.find("703_4530"),  # ~4.527 GeV
  DatasetManager.real_data.find("703_4575"),  # ~4.575 GeV
  DatasetManager.real_data.find("703_4600"),  # ~4.600 GeV
  DatasetManager.real_data.find("706_4610"),  # ~4.612 GeV
  DatasetManager.real_data.find("706_4620"),  # ~4.628 GeV
  DatasetManager.real_data.find("706_4640"),  # ~4.641 GeV
  DatasetManager.real_data.find("706_4660"),  # ~4.661 GeV
  DatasetManager.real_data.find("706_4680"),  # ~4.682 GeV
  DatasetManager.real_data.find("706_4700"),  # ~4.699 GeV
  DatasetManager.real_data.find("707_4740"),  # ~4.740 GeV
  DatasetManager.real_data.find("707_4750"),  # ~4.750 GeV
  DatasetManager.real_data.find("707_4780"),  # ~4.781 GeV
  DatasetManager.real_data.find("707_4840"),  # ~4.843 GeV
  DatasetManager.real_data.find("707_4914"),  # ~4.918 GeV
  DatasetManager.real_data.find("707_4946")   # ~4.951 GeV
]

# Matched inclusive MC samples (available subsets of the scan points)
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946")
]

# Decay card for the representative signal mode (EvtGen format):
#   psi(4260) -> omega chi_c1, omega -> pi+ pi- pi0, chi_c1 -> gamma J/psi,
#   J/psi -> e+ e-   (the same final state / selection also covers chi_c2 and the mu+mu- channel)
decay_card_omega_chic1 = <<~DECAYCARD
    Decay psi(4260)
    1.000 omega chi_c1 PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC generated at every energy point for the representative omega chi_c1 mode
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omega_chic1_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_omega_chic1
  config.cross_section = :default
end

### Event selection (BOSS) ###
# A single Algorithm (and a single selection) serves both chi_c1 / chi_c2 and both J/psi
# lepton channels, since they share the identical final state and selection.
alg_name = "OmegaChicJ"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})   # nominal CMS energy (per-point energy handled by the framework)
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:helix_correction, "Helix-parameter correction applied to all charged tracks before the 5C kinematic fit; efficiency difference with/without the correction estimated by re-running the BOSS selection.")

event_selection = Selection.new
event_selection.select_track {                    # Charged track selection
                    cos_theta  0.93               # |cos(theta)| < 0.93
                    Vz         10.0               # |Vz| < 10 cm
                    Vr         10.0               # Vr < 1 cm (= 10 mm) in transverse plane
                    nChrp      "==2"              # exactly two positive tracks
                    nChrn      "==2"              # exactly two negative tracks
                    nNet       "==0"              # net charge zero
                }
                .select_photon {                  # Photon selection
                    tdc_emc_start     0
                    tdc_emc_end       14
                    energyThreshold_b 0.025       # 25 MeV in the barrel
                    energyThreshold_e 0.050       # 50 MeV in the endcap
                    angle_to_track    10.0        # > 10 degrees from any charged track
                    nGam              ">=3"       # at least three photons
                }
                .pid(method: :probability) {      # Particle identification
                    # tracks with p > 1.0 GeV treated as leptons; electron if EMC energy > 1.0 GeV, else muon
                    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                   treat_as_electron_if_energy_above: 1.0
                    nlp ">=1"                     # at least one positive lepton
                    nlm ">=1"                     # at least one negative lepton
                }
                .remove([:lp <= :chrgp, :lm <= :chrgn])  # take the leptons out of the charged lists
                .assign({:chrgp => :pip, :chrgn => :pim})  # remaining tracks assigned to pi+ / pi-
                .kalman_kinematic_fit([:gamma, :gamma]) {  # reconstruct pi0 from two photons (1C)
                    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                    chi2_cut 25                   # chi2 < 25
                    npi0 ">=1"                    # at least one pi0 candidate
                }
                # Final 5C kinematic fit: 4C four-momentum conservation + the mass-constrained pi0
                # (tight published chi2 < 60 and the omega [0.75,0.81] / J/psi [3.08,3.12] GeV/c^2
                #  windows are applied later at ROOT level)
                .kinematic_fit([:pip, :pim, :lp, :lm, :gamma, :pi0]) {
                    nominal                       # nominal fit -> corrected four-momenta are stored
                    constrain_four_momentum       # 4C energy-momentum constraint
                    chi2_cut 200                  # loose chi2 cut in BOSS
                }

my_algorithm.with_decay_card(decay_card_omega_chic1).apply(event_selection)

# Execute on real data, matched inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs)