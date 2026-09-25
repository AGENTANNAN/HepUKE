# =============================================================================
#  J/psi event-number measurement by inclusive J/psi counting
#    1) main selection : count inclusive (hadronic) J/psi events @ 3.097 GeV
#    2) efficiency     : soft-pion tag of psi(2S) -> pi+ pi- J/psi @ 3.686 GeV
#  No particle identification and no kinematic fit are used in either selection.
# =============================================================================

### -------------------------------- Datasets --------------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data, 225.3 M events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # corresponding J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(2S) (3.686 GeV) real data, 106 M events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding psi(2S) inclusive MC

### ------------------------------- Decay cards ------------------------------- ###
# psi(2S) -> pi+ pi- J/psi, with J/psi -> anything.
# The J/psi is left undecayed in the card so that the generator's default
# (inclusive) decay table governs its decay.
decay_card_psip_pipi_jpsi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay
  End
DECAYCARD

# Direct J/psi -> anything : the J/psi is the top mother (produced directly) and
# is left to decay inclusively through the generator's default decay table.
decay_card_jpsi_direct = <<~DECAYCARD
  # J/psi -> anything (inclusive decay via the default decay table)
  End
DECAYCARD

### ------------------------------ Exclusive MC ------------------------------- ###
exMC_psip_pipi_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_to_pipi_jpsi_incl"  # 1 M-event signal MC
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psip_pipi_jpsi
  config.cross_section   = :default
end

exMC_jpsi_direct = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_direct_incl"        # 1 M-event direct J/psi MC
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_jpsi_direct
  config.cross_section   = :default
end

### ================== Algorithm 1 : inclusive J/psi counting ================ ###
alg_name = "JpsiInclusiveCount"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})            # J/psi CMS energy (GeV)
   .set_alias({"std::vector<double>" => "Vdouble"})

count_selection = Selection.new
count_selection.select_track {          # good charged tracks
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        15.0                    # |Vz| < 15 cm
      Vr        1.0                     # Vr < 1 cm
      momentum  2.0                     # p < 2.0 GeV/c
      nTot      ">=2"                   # at least two good charged tracks
    }
   .select_photon {                     # good photons (no multiplicity requirement)
      tdc_emc_start     0               # EMC timing 0 < T
      tdc_emc_end       15              # EMC timing T < 15
      energyThreshold_b 0.025           # > 25 MeV in the barrel
      energyThreshold_e 0.050           # > 50 MeV in the endcap
    }

# Event-level criteria that have no dedicated DSL primitive are captured as notes
# for the systematic-uncertainty stage.
alg.note(:visible_energy_cut,
         "event-level requirement E_vis > 1.0 GeV, where E_vis is the total visible
          energy summed over all good charged tracks and photons; applied after the
          track and photon selection")
   .note(:bhabha_dimuon_veto,
         "for events containing exactly two charged tracks, require both track
          momenta p < 1.5 GeV/c AND each track's EMC deposit < 1.0 GeV, to reject
          e+e- -> e+e- (Bhabha) and e+e- -> mu+mu- backgrounds; this condition is
          conditional on the two-track topology and cannot be written as a flat
          selection cut")

alg.with_decay_card(decay_card_jpsi_direct).apply(count_selection)

### ============== Algorithm 2 : psi(2S) -> pi+ pi- J/psi soft tag =========== ###
alg_name_eff = "PsipToPipiJpsiTag"
alg_eff = Algorithm.new(alg_name_eff)
alg_eff.set_header(["#{alg_name_eff}Alg/#{alg_name_eff}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})       # psi(2S) CMS energy (GeV)

eff_selection = Selection.new
eff_selection.select_track {            # soft-pion selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        15.0                    # |Vz| < 15 cm
      Vr        1.0                     # Vr < 1 cm
      momentum  0.4                     # p < 0.4 GeV/c (soft pions)
      nChrp     "==1"                   # exactly one pi+
      nChrn     "==1"                   # exactly one pi-
      nNet      "==0"                   # net charge zero
    }

# Detection-efficiency determination and the luminosity measurement are global,
# data-driven procedures without a dedicated DSL primitive.
alg_eff.note(:recoil_mass_fit,
             "detection efficiency is extracted from the psi(2S) -> pi+ pi- J/psi
              sample by fitting the J/psi recoil mass against the pi+ pi- pair;
              the result is corrected by the psi(2S)/J/psi inclusive-MC ratio
              (ratio of psi(2S) to J/psi inclusive-MC events passing the selection)")
       .note(:luminosity_from_gg,
             "integrated luminosity is measured independently from the QED process
              e+e- -> gamma gamma")

alg_eff.with_decay_card(decay_card_psip_pipi_jpsi).apply(eff_selection)

### -------------------------------- Execution -------------------------------- ###
root_files_count = alg.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_direct])
root_files_eff   = alg_eff.execute_on([psip_data, psip_incMC, exMC_psip_pipi_jpsi])