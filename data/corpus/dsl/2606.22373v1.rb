### dataset description ###
# Main dataset: J/psi at sqrt(s) = 3.097 GeV, 10.087 x 10^9 J/psi events
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Signal decay card: J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> e+- mu-+
decay_card_signal = <<~DECAYCARD
    Alias eta_prime  eta'

    Decay J/psi
    1.000  gamma  eta_prime          HELAMP 1 0 -1 0;
    Enddecay

    Decay eta_prime
    1.000  pi+  pi-  eta             PHSP;
    Enddecay

    Decay eta
    0.500  e+   mu-                  PHSP;
    0.500  e-   mu+                  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_jpsi_gam_etap_pipieta_emu"
  config.related_dataset = jpsi_data
  config.events = 500000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### event selection (BOSS) ###
alg_name = "EtaToEMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93        # |cos(theta)| < 0.93
                  Vz          10.0        # |Vz| < 10 cm
                  Vr          1.0         # |Vxy| < 1 cm
                  nChrp       "==2"       # two positive tracks (pi+, and e+ or mu+)
                  nChrn       "==2"       # two negative tracks (pi-, and e- or mu-)
                  nNet        "==0"       # net charge = 0
                }
               .select_photon {
                  tdc_emc_start    0      # EMC timing lower edge (0 ns)
                  tdc_emc_end      14     # EMC timing upper edge (~700 ns)
                  angle_to_track   10.0   # min angle to nearest charged track (deg)
                  energyThreshold_b 0.025 # barrel photon energy > 25 MeV
                  energyThreshold_e 0.050 # end-cap photon energy > 50 MeV
                  nGam             ">=1"  # at least one photon
                }
               .pid(method: :chi2_sum) {
                  # Combined chi2 PID (dE/dx + TOF + EMC) for pi+, pi-, e+/-, mu+/-,
                  # optimised together with the 4C kinematic fit.
                  identify :pion
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  npip   "==1"
                  npim   "==1"
                  nlp    "==1"
                  nlm    "==1"
               }
               # Nominal 4C kinematic fit for J/psi -> gamma pi+ pi- e+- mu-+
               .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 36            # min chi2_4C+PID < 36 (Punzi-optimised)
               }
               # Competing hypothesis: J/psi -> gamma pi+ pi- e+ e-
               # Stores chi2 for later cut chi2_4C+PID(pipiemu) < chi2_4C+PID(pipiee).
               .kinematic_fit([:gamma, :pip, :pim, :lp, :lp]) {
                  constrain_four_momentum
               }
               # Competing hypothesis: J/psi -> gamma pi+ pi- pi+ pi-
               # Stores chi2 for later cut chi2_4C+PID(pipiemu) < chi2_4C+PID(pipipipi).
               .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
                  constrain_four_momentum
               }

my_algorithm
  .note(:combined_chi2_PID, "PID chi2 (dE/dx + TOF + EMC) added to the 4C kinematic-fit chi2; combination minimising chi2_4C+chi2_PID(pi+ pi- e+/- mu-/+) is retained, chi2_4C+PID(pi+ pi- e+/- mu-/+) < 36")
  .note(:hypothesis_veto_ee, "veto: chi2_4C+PID(pi+ pi- e+/- mu-/+) < chi2_4C+PID(pi+ pi- e+ e-) required to suppress radiative eta' -> pi+ pi- e+ e- background")
  .note(:hypothesis_veto_4pi, "veto: chi2_4C+PID(pi+ pi- e+/- mu-/+) < chi2_4C+PID(pi+ pi- pi+ pi-) required to suppress eta' -> pi+ pi- pi+ pi- background")
  .note(:electron_EoverP, "additional electron identification: EMC energy deposit over track momentum E/p > 0.8 required for the e+/- candidate")
  .note(:muon_EMC_deposit, "additional muon identification: EMC energy deposit of the mu-/+ candidate required to lie in (0.1, 0.3) GeV")
  .note(:etap_mass_window, "M(pi+ pi- e+/- mu-/+) restricted to the eta' signal region (0.946, 0.970) GeV/c^2 (~ 3 sigma)")
  .note(:eta_signal_window, "M(e+/- mu-/+) signal window (0.538, 0.558) GeV/c^2 (~ 3 sigma of the eta mass resolution)")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
