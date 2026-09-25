### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")        # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # Corresponding inclusive MC sample

# Decay card for the signal process ψ(3686) → γ χ_c0 → γ π+π−π+π− (EvtGen format, PHSP generator)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process: 1,000,000 events, PHSP
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3686_gamma_chic0_4pi"
  config.related_dataset = psip_data     # Associated real dataset
  config.events = 1000000                # 1,000,000 events
  config.decay_card = decay_card_signal  # Decay card defining the γ 4π process
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "GammaChic0To4Pi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})  # ECMS = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain (called until the end of the final kinematic fit)
event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
                  cos_theta  0.93              # |cosθ| < 0.93
                  Vz         10.0              # |Vz| < 10 cm
                  Vr         1.0               # Vr < 1 cm
                  nChrp      "==2"             # exactly two positive tracks
                  nChrn      "==2"             # exactly two negative tracks
                  nNet       "==0"             # net charge zero
                }
               .select_photon {                # Photon selection
                  tdc_emc_start  0             # TDC start time
                  tdc_emc_end    14            # TDC end time
                  angle_to_track 10.0          # angle to nearest charged track > 10 degrees
                  energyThreshold_b 0.025      # 25 MeV barrel energy threshold
                  energyThreshold_e 0.050      # 50 MeV endcap energy threshold
                  nGam  ">=1"                  # at least one photon
                }
               .pid(method: :probability) {    # PID by the probability method
                  prob_cut 0.001               # PID probability > 0.001
                  identify :pion, against: [:kaon, :proton]  # π+ and π− vs K and p
                  npip  "==2"                  # two π+ identified
                  npim  "==2"                  # two π− identified
                }
                # Nominal 4C kinematic fit to γ π+π−π+π− with a common vertex constraint on the four pions
               .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
                  nominal                      # nominal fit — corrected four-momenta are used
                  vertex_fit([1, 2, 3, 4])     # common vertex fit on the four pions (indices 1–4)
                  constrain_four_momentum      # 4C energy-momentum constraint
                  chi2_cut 40                  # χ² < 40
                }
                # Alternative 3C fit (omitting photon energy), stored for signal extraction
               .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
                  constrain_three_momentum     # 3C-type constraint, chi2 stored for later use
                }

# Capture BOSS-side procedures / criteria that the DSL cannot express
my_algorithm
  .note(:background_veto, "all π+π− recoil masses (M_recoil² = (P_CMS − p_ππ)²) required to lie " \
        "outside the J/ψ region 3.0–3.2 GeV/c², applied combinatorially over all π+π− pairs, " \
        "to suppress ψ(3686) → π+π−J/ψ background")
  .note(:background_veto, "all γπ+π− recoil masses required to lie outside the J/ψ region " \
        "3.0–3.2 GeV/c², applied over all γπ+π− combinations")
  .note(:background_veto, "cosθ of every π+π− combination required within [−0.999, 0.988] " \
        "to reject γ → e+e− conversion background")
  .note(:fsr_correction, "FSR correction factor f_FSR = 2.00 ± 0.02 applied to MC; the ±0.02 " \
        "uncertainty is treated as a systematic source")
  .note(:helix_correction, "helix-parameter correction applied to charged tracks in MC before the " \
        "4C kinematic fit")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])