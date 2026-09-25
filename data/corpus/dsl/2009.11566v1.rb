# Trigger efficiency study with 2018 J/psi data — detector performance paper
# This analysis measures trigger efficiencies using control samples, not decay branching fractions.

psip_data = DatasetManager.real_data.find("708_3097")
psip_incMC = DatasetManager.inclusive_mc.find("708_3097")

algorithm = Algorithm.new("TriggerEffStudy")
algorithm.set_header(["TriggerEffStudyAlg/TriggerEffStudy.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .note(:scope,
            "Detector trigger-efficiency study, not a physics decay analysis. " \
            "Control samples selected: Bhabha (e+e-), dimuon (mu+mu-), " \
            "and generic hadronic events from J/psi data.")
          .note(:bhabha_selection,
            "Bhabha: 2 EMC clusters with opening > 166deg, " \
            "|E_emc(e+)+E_emc(e-)-3.097|/3.097 <= 0.10. " \
            "2 opp-charge good tracks, opening > 175deg.")
          .note(:dimuon_selection,
            "Dimuon: 2 opp-charge good tracks, opening >= 178deg, " \
            "each track p < 2 GeV/c, E_EMC < 0.7 GeV. " \
            "Total 4-momentum (E,px,py,pz) in [(2.8,3.3),(-0.1,0.1),(-0.1,0.1),(-0.2,0.2)].")
          .note(:hadronic_selection,
            "Hadronic: >=2 good tracks in MDC. If exactly 2 tracks, " \
            "opening < 170deg to suppress Bhabha/dimuon.")
          .note(:trigger_efficiency_formula,
            "Efficiency = N(sel & trig_cond/chan) / N(sel). " \
            "Global efficiency from inclusion-exclusion over 3 groups of trigger channels. " \
            "Clopper-Pearson 68.27% CL intervals.")

# Minimal shared selection: good tracks with polar angle cut
event_selection = Selection.new
event_selection.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
  nNet  "==0"
}

algorithm.apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC])