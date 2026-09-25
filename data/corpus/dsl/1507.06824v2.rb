# Measurement of Azimuthal Asymmetries in Inclusive Charged Dipion Production
# e+ e- -> pi pi X at sqrt(s) = 3.65 GeV (Collins fragmentation function)

### Dataset description ###
# 62 pb^-1 continuum data collected at 3.65 GeV.
data_3650  = DatasetManager.real_data.find("703_3650")
incMC_3650 = DatasetManager.inclusive_mc.find("703_3650")

# Inclusive continuum e+e- -> qqbar (u,d,s); no exclusive final state, no top-mother
# resonance. Use the ConExc/psi(4260) convention just to keep the DSL template happy.
decay_card_qq = <<~DECAYCARD
  Decay psi(4260)
  1.000 u anti-u                            PHSP;
  Enddecay

  End
DECAYCARD

exMC_qq = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "collins_uds_continuum"
  config.related_dataset = data_3650
  config.events         = 500000
  config.decay_card     = decay_card_qq
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "CollinsInclDipion"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.650]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nTot      ">=3"
    }
    .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon, :proton]
      identify :kaon, against: [:pion, :proton]
      npip ">=1"
      npim ">=1"
    }

alg.note(:electron_veto,
         "Reject events containing any identified electron; electron ID requires " \
         "L(e) > 0.001 and L(e)/(L(e)+L(pi)+L(K)) > 0.8.")
   .note(:visible_energy_cut,
         "Total visible energy (all reconstructed charged tracks + photons) > 1.5 GeV " \
         "to suppress QED tau+tau- and beam-gas backgrounds.")
   .note(:emc_timing,
         "EMC cluster timing delay from event start time <= 700 ns.")
   .note(:pi_pid,
         "Charged pion ID: L(pi) > 0.001, L(pi) > L(K), L(pi) > L(p).")
   .note(:kaon_pid,
         "Charged kaon ID (used to estimate mis-ID contamination): L(K) > 0.001, " \
         "L(K) > L(pi), L(K) > L(p).")
   .note(:dipion_z_range,
         "Fractional energy z_i = 2 E_pi / sqrt(s) required to be in [0.2, 0.9] for both pions.")
   .note(:opening_angle,
         "Opening angle between the two charged pion candidates > 120 degrees to select " \
         "back-to-back pairs (proxy for the qqbar axis in absence of clear jet structure at 3.65 GeV).")
   .note(:pair_combinatorics,
         "Two pions per pair are labelled randomly h1 and h2; if more than two pions are " \
         "present in an event, all combinations are formed and each pion may enter multiple pairs.")

alg.with_decay_card(decay_card_qq).apply(sel)
alg.execute_on([data_3650, incMC_3650, exMC_qq])
