# frozen_string_literal: true
### Dataset description ###
# J/psi (3.097 GeV) detector trigger-efficiency study: 2018 real data + matching inclusive MC.
# The 2018 runs belong to round 11/12 (runs 52940-54976, 55861-56546, 56788-59015).
jpsi_data  = DatasetManager.real_data.find("708_3097")     # real J/psi data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # matching inclusive MC sample

### Decay cards ###
# This is a control-sample trigger-efficiency study rather than a decay analysis, so no
# exclusive MC is generated. One decay card per algorithm is still required to seed the
# kinematic-variable header; each card reflects the control-sample topology.
decay_card_bhabha = <<~DECAYCARD
  Decay psi(4260)
  1.0000 e+ e- PHSP;
  Enddecay
  End
DECAYCARD

decay_card_dimuon = <<~DECAYCARD
  Decay psi(4260)
  1.0000 mu+ mu- PHSP;
  Enddecay
  End
DECAYCARD

# Representative generic hadronic final state (placeholder used only for the header).
decay_card_hadronic = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

### Common preselection ###
# good charged tracks: |cos(theta)| < 0.93, |Vz| < 10 cm, Vr < 1 cm,
# at least one positive and one negative track, net charge zero.
# It is written inside each control-sample select_track block because the required
# track multiplicity differs between the three control samples.

### --- Bhabha control sample: e+e- -> e+e- --- ###
alg_bhabha = Algorithm.new("TrgEffBhabha")
alg_bhabha.set_header(["TrgEffBhabhaAlg/TrgEffBhabha.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .note(:cluster_opening_angle, "the two EMC clusters of the Bhabha candidate must have an opening angle > 166 deg; inter-cluster opening-angle cuts are not expressible in the DSL and are applied in ROOT")
          .note(:cluster_energy_sum, "the Bhabha candidate must satisfy |E_emc(e+)+E_emc(e-)-3.097|/3.097 <= 0.10; summing EMC cluster energies is not expressible in the DSL and is evaluated in ROOT")
          .note(:track_opening_angle, "the two opposite-charge good tracks must have an opening angle > 175 deg; inter-track opening-angle cuts are not expressible in the DSL and are applied in ROOT")
          .note(:trigger_condition, "the per-event Bhabha-channel trigger decision must be evaluated and stored; the efficiency N(trigger & selected)/N(selected) is formed in the ROOT analysis")

bhabha_selection = Selection.new
  .select_track {
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz 10.0            # |Vz| < 10 cm (beam direction)
    Vr 1.0             # Vr < 1 cm (transverse plane)
    nChrp "==1"        # exactly one positive track
    nChrn "==1"        # exactly one negative track
    nNet  "==0"        # net charge zero
  }
  .select_photon {
    nGam "==2"         # exactly two EMC clusters (angle / energy criteria handled in ROOT)
  }

alg_bhabha.with_decay_card(decay_card_bhabha).apply(bhabha_selection)
root_files_bhabha = alg_bhabha.execute_on([jpsi_data, jpsi_incMC])

### --- Dimuon control sample: e+e- -> mu+mu- --- ###
alg_dimuon = Algorithm.new("TrgEffDimuon")
alg_dimuon.set_header(["TrgEffDimuonAlg/TrgEffDimuon.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .note(:track_opening_angle, "the two opposite-charge good tracks must have an opening angle >= 178 deg; inter-track opening-angle cuts are not expressible in the DSL and are applied in ROOT")
          .note(:track_momentum_emc, "each of the two tracks must satisfy p < 2 GeV/c and E_EMC < 0.7 GeV; per-track momentum and track-associated EMC energy cuts are not expressible in the DSL and are applied in ROOT")
          .note(:four_momentum_window, "the total four-momentum must satisfy (E,px,py,pz) in (2.8,3.3),(-0.1,0.1),(-0.1,0.1),(-0.2,0.2) GeV; summing the pair four-momentum is not expressible in the DSL and is applied in ROOT")
          .note(:trigger_condition, "the per-event dimuon-channel trigger decision must be evaluated and stored; the efficiency N(trigger & selected)/N(selected) is formed in the ROOT analysis")

dimuon_selection = Selection.new
  .select_track {
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz 10.0            # |Vz| < 10 cm
    Vr 1.0             # Vr < 1 cm
    nChrp "==1"        # exactly one positive track
    nChrn "==1"        # exactly one negative track
    nNet  "==0"        # net charge zero
  }

alg_dimuon.with_decay_card(decay_card_dimuon).apply(dimuon_selection)
root_files_dimuon = alg_dimuon.execute_on([jpsi_data, jpsi_incMC])

### --- Generic hadronic control sample --- ###
alg_hadronic = Algorithm.new("TrgEffHadronic")
alg_hadronic.set_header(["TrgEffHadronicAlg/TrgEffHadronic.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .note(:two_track_opening_angle, "for events with exactly two good tracks the opening angle is required to be < 170 deg to suppress Bhabha and dimuon contamination; this two-track opening-angle veto is not expressible in the DSL and is applied in ROOT")
            .note(:trigger_condition, "the per-event hadronic-channel trigger decision must be evaluated and stored; the efficiency N(trigger & selected)/N(selected) is formed in the ROOT analysis")

hadronic_selection = Selection.new
  .select_track {
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz 10.0            # |Vz| < 10 cm
    Vr 1.0             # Vr < 1 cm
    nChrp ">=1"        # at least two good MDC tracks in total ...
    nChrn ">=1"        # ... (>=1 of each charge, with net charge zero)
    nNet  "==0"
  }

alg_hadronic.with_decay_card(decay_card_hadronic).apply(hadronic_selection)
root_files_hadronic = alg_hadronic.execute_on([jpsi_data, jpsi_incMC])