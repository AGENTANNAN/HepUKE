### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding inclusive MC sample

# Decay card: psi(2S) -> anti-p p+ phi (phase space), phi -> K+ K- (VSS)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000  anti-p-  p+  phi    PHSP;
    Enddecay

    Decay phi
    1.0000  K+  K-    VSS;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_ppbarphi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name     = "PsipToPpbarPhi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})       # CMS energy = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {              # Charged-track selection
                  cos_theta  0.93           # |cos(theta)| < 0.93
                  Vz         100.0          # |Vz| < 100 cm
                  Vr         10.0           # Vr < 10 (transverse plane)
                  nChrp      ">=2"          # at least two positive tracks
                  nChrn      ">=2"          # at least two negative tracks
                  nNet       "==0"          # net charge zero
                }
               .pid(method: :probability) {  # PID by the probability method
                  prob_cut   0.001           # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]   # p+ / p-bar vs K and pi
                  identify :kaon,   against: [:pion, :proton] # K+ / K- vs pi and p
                  nprp   ">=1"               # at least one proton
                  nprm   ">=1"               # at least one anti-proton
                  nkp    ">=1 || nkm >= 1"   # at least one kaon (K+ or K-)
                }
               .remove([:prp <= :chrgp, :prm <= :chrgn,   # drop identified (anti-)protons
                        :kp  <= :chrgp, :km  <= :chrgn])  # drop identified K+/K-
               .assign({:chrgp => :pip, :chrgn => :pim})  # remaining tracks treated as pions
                # 1C kinematic fit: psi(2S) -> p p-bar K+ K-, one kaon allowed to be missing.
                # When both kaons are detected, the combination with the smallest chi2 is chosen.
               .kinematic_fit([:prp, :prm, :kp, :km]) {
                  nominal
                  miss_track_of :km          # one kaon (K-) may be undetected
                  constrain_four_momentum    # constrain four-momentum to the CMS energy
                  chi2_cut 200               # loose chi2 cut (tight cut applied later in ROOT)
                }

# BOSS-side procedures with no direct DSL construct are preserved as notes.
my_algorithm
  .note(:background_veto, "events with M(p pi-) or M(p-bar pi+) within 3 MeV/c^2 of the " \
        "Lambda mass (1.1157 GeV/c^2, i.e. the [1.1127, 1.1187] GeV/c^2 window) are vetoed " \
        "to suppress Lambda -> p pi- and Lambda-bar -> p-bar pi+ contamination")
  .note(:fit_hypothesis, "the 1C kinematic fit admits either kaon (K+ or K-) as the missing " \
        "particle; the DSL emits a single nominal fit using K- as the missing track")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])