### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # inclusive MC for J/psi

# Decay card: J/psi -> K- Sigma0 anti-Xi+
# c.c. mode is included via BesEvtGen conjugation as usual.
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 K- Sigma0 anti-Xi+                    PHSP;
    Enddecay

    Decay anti-Xi+
    1.000 anti-Lambda0 pi+                      PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+                           HypWK;
    Enddecay

    Decay Sigma0
    1.000 Lambda0 gamma                         PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi-                                HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_jpsi_km_sigma0_xibarp"
    config.related_dataset = jpsi_data
    config.events         = 200000
    config.decay_card     = decay_card_signal
    config.cross_section  = :default
end

### Event selection ###
alg_name = "JpsiKmSigma0XibarP"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz        20.0    # |Vz| < 20 cm  (Lambda, Xi are long-lived)
                  Vr        100.0   # Vxy < 10 cm
                  nChrp     ">=2"   # anti-p and two pi+ (from Xi_bar+ -> Lambda_bar pi+, Lambda_bar -> pbar pi+)
                  nChrn     ">=1"   # K-
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # p and anti-p
                  identify :kaon,   against: [:proton, :pion] # K+ / K-
                  identify :pion,   against: [:kaon, :proton]
                  nprm ">=1"        # one anti-proton
                  nkm  ">=1"        # one K-
                  npip ">=2"        # two pi+
                }
               # K- originates from IP directly: additional tight vertex requirements
               .for_each(:km) {
                  where { Vxy > 1.0 }         # |Vxy| < 1 cm
                  remove
               }
               .for_each(:km) {
                  where { abs(Vz) > 5.0 }     # |Vz| < 5 cm
                  remove
               }
               # Reconstruct anti-Lambda from (anti-p, pi+) via secondary vertex fit
               .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # Reconstruct anti-Xi+ from (anti-Lambda, pi+) via secondary vertex fit
               .secondary_vertex_fit([:Lambda_bar, :pip]) {
                  build_virtual_particle(:anti_Xi_p).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # Partial reconstruction: tag K- and anti-Xi+; Sigma0 is the missing recoil particle
               .partial_miss([2]) {           # Sigma0 (rec_id 2) treated as missing
                  require_recoil_mass 1.174, 1.204   # |M_recoil(K- anti-Xi+) - m_Sigma0| < 15 MeV/c^2
               }

alg.note(:xibar_selection,
         "anti-Xi+ mass window: |M - m_Xi| < 10 MeV/c^2 and decay length > 0. " \
         "Lambda/Xi peaking suppression: |(M(Lbar pi+) - M(Xibar+)) - (M(pbar pi+) - M(Lbar))| < 4 MeV/c^2.")
   .note(:recoil_kaon_cut,
         "M_recoil(K-) < 2.58 GeV/c^2 to discard events with low-momentum kaons (poor resolution).")
   .note(:sigma0_mass_constraint,
         "1-C kinematic fit constraining missing mass to nominal Sigma0 mass, applied to " \
         "improve mass resolution only (no chi2 cut).")

alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
