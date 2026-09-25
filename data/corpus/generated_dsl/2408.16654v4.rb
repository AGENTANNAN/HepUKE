# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Corresponding inclusive MC sample

# Decay card for the signal process J/psi -> Xi0 anti-Xi0,
# Xi0 -> Lambda pi0, anti-Xi0 -> anti-Lambda gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi0 anti-Xi0          PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0           PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 gamma    PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+           HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_xi0_to_lambdapi0_xibar0_to_lambdabargamma"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Xi0ToLambdaGamma"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                       # Charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        100.0                   # |Vz| < 100 mm
      Vr        10.0                    # Vr < 10 mm
      nChrp     ">=2"                   # At least two positive tracks
      nChrn     ">=2"                   # At least two negative tracks
      nNet      "==0"                   # Net charge zero
  }
  .select_photon {                      # Photon selection
      tdc_emc_start     0               # EMC TDC start
      tdc_emc_end       14              # EMC TDC end
      angle_to_track    10.0            # Angle to nearest charged track > 10 deg
      energyThreshold_b 0.025           # Barrel energy > 25 MeV
      energyThreshold_e 0.050           # Endcap energy > 50 MeV
      nGam              ">=2"           # At least two photons
  }
  .pid(method: :probability) {          # Particle identification (probability method)
      prob_cut 0.001                    # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]  # identify p and p-bar vs K and pi
      nprp     ">=1"                    # At least one proton
      nprm     ">=1"                    # At least one anti-proton
  }
  .remove([:prp <= :chrgp])             # Remove identified protons from positive list
  .remove([:prm <= :chrgn])             # Remove identified anti-protons from negative list
  .assign({:chrgp => :pip, :chrgn => :pim})  # Remaining tracks assigned as pi+ / pi-
  .secondary_vertex_fit([:prp, :pim]) { # Lambda -> p pi- from a common secondary vertex
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) { # anti-Lambda -> p-bar pi+ from a common secondary vertex
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {  # 1C Kalman fit to reconstruct pi0 -> gamma gamma
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25                       # chi2 < 25
      npi0     ">=1"                    # At least one pi0
  }
  .kinematic_fit([:Lambda, :pi0, :Lambda_bar, :gamma]) {  # Nominal 4C kinematic fit to Lambda pi0 anti-Lambda gamma
      nominal                           # Mark as the nominal fit
      constrain_four_momentum           # 4-momentum conservation (pi0 mass fixed by Kalman step)
      chi2_cut 40                       # chi2 < 40
  }
  .kinematic_fit([:Lambda, :pi0, :Lambda_bar, :gamma]) {  # Competing Sigma0-hypothesis fit
      constrain_four_momentum           # same constraint, no nominal flag (chi2 stored for ROOT-level veto)
  }

# Attach the decay card and render the selection into the BOSS algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])