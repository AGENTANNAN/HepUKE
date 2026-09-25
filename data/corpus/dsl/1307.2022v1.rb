### Dataset description ###
# psi(3770) data at sqrt(s) = 3.773 GeV
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Continuum data at sqrt(s) = 3.650 GeV
cont3650_data  = DatasetManager.real_data.find("712_3650")
cont3650_incMC = DatasetManager.inclusive_mc.find("712_3650")

# Decay card for Bhabha scattering signal MC: e+ e- -> (gamma) e+ e-
# Use the Babayaga generator (BESIII convention token). The Bhabha final state
# is directly produced from the e+ e- initial system, so we adopt the KKMC top
# mother convention with psi(4260).
decay_card_bhabha = <<~DECAYCARD
    Decay psi(4260)
    1.0000  e+ e-                              BABAYAGA;
    Enddecay

    End
DECAYCARD

# Exclusive Bhabha MC generated at both energy points (3.650 GeV and 3.773 GeV)
exMC_bhabha_set = DatasetManager.create_exclusive_mc_for([cont3650_data, psi3770_data]) do |config|
  config.sample_name    = "bhabha_lumi_exclusive_mc"
  config.events         = 400_000
  config.decay_card     = decay_card_bhabha
  config.cross_section  = :default
end
exMC_bhabha_set.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) — Bhabha scattering event selection for luminosity ###

# ---- Algorithm for the psi(3770) data set (sqrt(s) = 3.773 GeV) ----
alg_3773 = Algorithm.new("BhabhaLumi3773")
alg_3773.set_header(["BhabhaLumi3773Alg/BhabhaLumi3773.h"])
        .set_constant({"ECMS" => [:double, 3.773]})

sel_3773 = Selection.new
sel_3773.select_track {
           cos_theta 0.80        # |cos(theta)| < 0.80 -> track hits EMC barrel
           Vr         10.0       # Rxy < 1 cm (10 mm) closest approach in xy
           Vz         5.0        # |Vz| < 5 cm at closest approach
           nChrp     "==1"       # exactly one positively charged track
           nChrn     "==1"       # exactly one negatively charged track
           nNet      "==0"       # total charge zero
         }
        .pid(method: :probability) {
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
         }
        .for_each(:chrgp) {
           define(:eEmcp) { emc_energy }
           where { eEmcp > 1.0 } # E_EMC of positron > 1.0 GeV (reject mu+mu- events)
         }
        .for_each(:chrgn) {
           define(:eEmcn) { emc_energy }
           where { eEmcn > 1.0 } # E_EMC of electron > 1.0 GeV (reject mu+mu- events)
         }

alg_3773.note(:cosmic_veto,
              "Reject cosmic rays: the momenta of the two charged tracks must " \
              "NOT both exceed E_beam + 0.15 GeV (E_beam = ECMS/2).")
        .note(:jpsi_background_veto,
              "Suppress (gamma)J/psi, (gamma)psi(3686) and psi(3770) -> (gamma)J/psi X " \
              "contamination by requiring |p1| + |p2| > 0.9 * ECMS.")
        .note(:signal_region,
              "Signal extraction uses delta_phi = |phi1 - phi2| - 180 deg between the " \
              "two EMC clusters; signal region and sideband defined on this variable " \
              "for background subtraction. Applied in ROOT.")
        .note(:trigger_efficiency,
              "Trigger efficiency for e+e- -> (gamma) e+ e- taken as 100% (stat. error < 0.1%).")
        .note(:generator, "Signal Bhabha MC generated with Babayaga, |cos theta| < 0.83.")

alg_3773.with_decay_card(decay_card_bhabha).apply(sel_3773)
alg_3773.execute_on([psi3770_data, psi3770_incMC, exMC_bhabha_set[1]])

# ---- Algorithm for the continuum data set (sqrt(s) = 3.650 GeV) ----
alg_3650 = Algorithm.new("BhabhaLumi3650")
alg_3650.set_header(["BhabhaLumi3650Alg/BhabhaLumi3650.h"])
        .set_constant({"ECMS" => [:double, 3.650]})

sel_3650 = Selection.new
sel_3650.select_track {
           cos_theta 0.80
           Vr         10.0
           Vz         5.0
           nChrp     "==1"
           nChrn     "==1"
           nNet      "==0"
         }
        .pid(method: :probability) {
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
         }
        .for_each(:chrgp) {
           define(:eEmcp) { emc_energy }
           where { eEmcp > 1.0 }
         }
        .for_each(:chrgn) {
           define(:eEmcn) { emc_energy }
           where { eEmcn > 1.0 }
         }

alg_3650.note(:cosmic_veto,
              "Reject cosmic rays: momenta of the two tracks must NOT both exceed " \
              "E_beam + 0.15 GeV (E_beam = ECMS/2).")
        .note(:jpsi_background_veto,
              "Suppress residual (gamma)J/psi / (gamma)psi(3686) backgrounds by " \
              "requiring |p1| + |p2| > 0.9 * ECMS.")
        .note(:signal_region,
              "Signal extraction via delta_phi = |phi1 - phi2| - 180 deg between " \
              "the two EMC clusters; sideband subtraction applied in ROOT.")
        .note(:trigger_efficiency,
              "Trigger efficiency for e+e- -> (gamma) e+ e- taken as 100% (stat. error < 0.1%).")
        .note(:generator, "Signal Bhabha MC generated with Babayaga, |cos theta| < 0.83.")

alg_3650.with_decay_card(decay_card_bhabha).apply(sel_3650)
alg_3650.execute_on([cont3650_data, cont3650_incMC, exMC_bhabha_set[0]])
