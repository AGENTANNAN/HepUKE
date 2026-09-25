### Dataset preparation ###
# e+e- -> K+K- Born cross-section R-scan at sqrt(s) = 2.00-3.08 GeV (BOSS 713)
# 22 R-scan energy points (sample name convention: [BOSS]_[CMS energy in MeV])
rscan_names = ["713_2000", "713_2050", "713_2100", "713_2125", "713_2150",
               "713_2175", "713_2200", "713_2232", "713_2309", "713_2386",
               "713_2396", "713_2500", "713_2644", "713_2646", "713_2700",
               "713_2800", "713_2900", "713_2950", "713_2981", "713_3000",
               "713_3020", "713_3080"]
rscan_data  = rscan_names.map { |n| DatasetManager.real_data.find(n) }     # real data, all 22 points
rscan_incMC = rscan_names.map { |n| DatasetManager.inclusive_mc.find(n) }  # inclusive MC, all 22 points

# ConExc decay card for the ISR Born-cross-section measurement (mode 45 = K+K-).
# The DSL auto-detects the literal `ConExc` token, switches to the no-KKMC template
# and injects `Particle vpho <ECMS> 0.0` per energy point -- so `Particle vpho` is omitted.
decay_card_kkmc = <<~DECAYCARD
    Decay vpho
    1.00000 ConExc 45;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC generated at each of the 22 energy points
exMC_signal = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "kkmc_born_exclusive_mc"  # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_kkmc
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "KKmcBorn"
kk_alg = Algorithm.new(alg_name)
kk_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({"ECMS" => [:double, 3.08]})  # per-point vpho energy injected by the ConExc template
      .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # exactly one K+ and one K- candidate
        cos_theta 0.93                  # |cos(theta)| < 0.93
        Vz        10.0                  # |Vz| < 10 cm
        Vr        1.0                   # Vr < 1 cm
        nChrp     "==1"                 # exactly one positive charged track
        nChrn     "==1"                 # exactly one negative charged track
        nNet      "==0"                 # net charge zero
      }
      .pid(method: :probability) {       # kaon identification, no lepton requirement
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        nkp "==1"
        nkm "==1"
      }
      .kinematic_fit([:kp, :km]) {       # 4-momentum-constrained fit on the K+K- pair
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

kk_alg
  .note(:background_veto, "Additional optimized cuts suppress the e+e- and (gamma)mu+mu- backgrounds: energy-dependent E/p < 0.7-0.8, cos(theta) < 0.8 for the K+ and > -0.8 for the K-, K+K- opening angle > 179 deg in the CM frame, |Delta TOF| < 3 ns (cosmic rejection), and a +/-3 sigma momentum window on the K- track. The E/p thresholds are energy-dependent and were tuned per scan point; these cuts have no dedicated DSL method.")
  .with_decay_card(decay_card_kkmc)
  .apply(event_selection)

# Execute on real data, inclusive MC and the per-point exclusive signal MC
root_files = kk_alg.execute_on(rscan_data + rscan_incMC + exMC_signal)