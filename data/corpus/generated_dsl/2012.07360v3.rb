# =============================================================================
#  e+e- -> eta' pi+ pi-  :  Born cross-section scan via the ConExc generator
#  Mode I : eta' -> eta  pi+ pi- , eta -> gamma gamma
#  Mode II: eta' -> gamma pi+ pi-
#  BOSS release 713, R-scan energy points 2.000 - 3.080 GeV
# =============================================================================

### ---------------------------------------------------------------------------
### Dataset preparation
### ---------------------------------------------------------------------------

# R-scan energy points (MeV) covered by the scan (BOSS 713)
rscan_ecms = [2000, 2050, 2100, 2150, 2175, 2200, 2232, 2309, 2386, 2396,
              2500, 2644, 2646, 2700, 2800, 2900, 2950, 2981, 3000, 3020, 3080]

# Real data and the matching inclusive MC at every energy point
rscan_data  = rscan_ecms.map { |e| DatasetManager.real_data.find("713_#{e}") }
rscan_incMC = rscan_ecms.map { |e| DatasetManager.inclusive_mc.find("713_#{e}") }

# --- Decay card, Mode I: vpho -> eta' pi+ pi-, eta' -> eta pi+ pi-, eta -> gamma gamma ---
# The ConExc model supplies ISR (up to second order) and the measured sigma0(sqrt(s)).
# "Particle vpho" is intentionally omitted: the DSL injects it per energy point of the scan.
decay_card_modeI = <<~DECAYCARD
    Decay vpho
    1.0000 eta' pi+ pi- ConExc vhdr;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card, Mode II: vpho -> eta' pi+ pi-, eta' -> gamma pi+ pi- ---
decay_card_modeII = <<~DECAYCARD
    Decay vpho
    1.0000 eta' pi+ pi- ConExc vhdr;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal MC generated at every scan point, one sample per point
exMC_modeI = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "rscan_etapPiPi_modeI"
  config.events        = 200000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "rscan_etapPiPi_modeII"
  config.events        = 200000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode I : eta' -> eta pi+ pi-, eta -> gamma gamma
### ---------------------------------------------------------------------------
alg_name_modeI = "EtapPiPiEtaModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.080]})   # constant placeholder for the scan
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:isr_vp_correction,
               "ISR and vacuum-polarisation correction factors are taken from the ConExc generator log; they are applied for the Born cross-section extraction downstream in ROOT")

sel_modeI = Selection.new
sel_modeI.select_track {                      # charged tracks
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     ">=2"                    # >= 2 positive tracks
           nChrn     ">=2"                    # >= 2 negative tracks
         }
         .select_photon {                     # photons
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0             # min angle to nearest charged track
           energyThreshold_b 0.025            # E > 25 MeV (barrel)
           energyThreshold_e 0.050            # E > 50 MeV (endcap)
           nGam              ">=4"            # Mode I: >= 4 photons
         }
         .pid(method: :probability) {         # PID : probability method
           prob_cut 0.001
           identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K and p
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {    # eta -> gamma gamma (mass-constrained)
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 25
           neta     ">=1"
         }
         .kinematic_fit([:pip, :pim, :pip, :pim, :eta]) {   # nominal 4C fit
           nominal
           constrain_four_momentum
           chi2_cut 100
         }
         .kinematic_fit([:pip, :pim, :pip, :pim]) {         # competing hypothesis 2(pi+ pi-)
           constrain_four_momentum                          # chi2 stored, vetoed in ROOT
         }
         .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) {  # competing hypothesis 2(gamma pi+ pi-)
           constrain_four_momentum                                   # chi2 stored, vetoed in ROOT
         }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode II : eta' -> gamma pi+ pi-
### ---------------------------------------------------------------------------
alg_name_modeII = "EtapPiPiGammaModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.080]})   # constant placeholder for the scan
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:isr_vp_correction,
                "ISR and vacuum-polarisation correction factors are taken from the ConExc generator log; they are applied for the Born cross-section extraction downstream in ROOT")

sel_modeII = Selection.new
sel_modeII.select_track {                     # charged tracks
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     ">=2"                   # >= 2 positive tracks
            nChrn     ">=2"                   # >= 2 negative tracks
          }
          .select_photon {                    # photons
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0            # min angle to nearest charged track
            energyThreshold_b 0.025           # E > 25 MeV (barrel)
            energyThreshold_e 0.050           # E > 50 MeV (endcap)
            nGam              ">=3"           # Mode II: >= 3 photons
          }
          .pid(method: :probability) {        # PID : probability method
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K and p
          }
          .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma, :gamma]) {  # nominal 4C fit
            nominal
            constrain_four_momentum
            chi2_cut 50
          }
          .kinematic_fit([:pip, :pim, :pip, :pim]) {         # competing hypothesis 2(pi+ pi-)
            constrain_four_momentum                          # chi2 stored, vetoed in ROOT
          }
          .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) {  # competing hypothesis 2(gamma pi+ pi-)
            constrain_four_momentum                                   # chi2 stored, vetoed in ROOT
          }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

### ---------------------------------------------------------------------------
### Execution on real data, inclusive MC and exclusive signal MC at every point
### ---------------------------------------------------------------------------
root_files_modeI  = alg_modeI.execute_on(rscan_data + rscan_incMC + exMC_modeI)
root_files_modeII = alg_modeII.execute_on(rscan_data + rscan_incMC + exMC_modeII)