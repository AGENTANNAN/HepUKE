# 2205.04809v2: Luminosities and energies of e+e- collision data
# taken between sqrt(s)=4.61 GeV and 4.95 GeV at BESIII
#
# NOTE: This is a detector-performance paper (CMS energy and luminosity
# measurement), NOT a physics event-selection analysis. The CMS energies
# are measured via Lambda_c+ Lambda_c- partial reconstruction
# (Lambda_c+ -> p K- pi+), and luminosities via large-angle Bhabha
# scattering e+e- -> (gamma) e+e-. A di-photon cross-check is also
# performed. These are three independent measurements with no common
# "signal selection" in the BESIII DSL sense. This spec captures the
# common track-quality and PID requirements and documents the three
# measurement paths via notes.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# --- Decay card for Lambda_c+ Lambda_c- production ---
decay_card_lclc = <<~DECAYCARD
  Decay psi(4260)
  1.0 Lambda_c+ Lambda_c- PHSP;
  Enddecay
  End
DECAYCARD

# --- Data: 12 energy points, BOSS 706 (4.61-4.70 GeV) and 707 (4.74-4.95 GeV) ---
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

all_data = [data_4610, data_4620, data_4640, data_4660, data_4680, data_4700,
            data_4740, data_4750, data_4780, data_4840, data_4914, data_4946]

# --- Common track-quality selection (applies to all three measurements) ---
event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    identify :proton, against: [:kaon, :pion]
  end

# --- Common algorithm skeleton (multi-energy scan) ---
algorithm = Algorithm.new("LcEnergyLumiMeasurement")
algorithm.set_header(["LcEnergyLumiMeasurementAlg/LcEnergyLumiMeasurement.h"])
         .set_constant({ "ECMS" => [:double, 4.640] })  # placeholder; overridden per run

algorithm.note(:cms_energy_method,
  "CMS energies measured via Lambda_c+ Lambda_c- partial reconstruction: " \
  "Lambda_c+ -> p K- pi+. The invariant mass of partially reconstructed " \
  "Lambda_c pairs determines ECM per run. NOT expressible as a BOSS-side " \
  "kinematic fit selection — performed offline via ROOT analysis.")

algorithm.note(:luminosity_method,
  "Luminosity measured via large-angle Bhabha scattering e+e- -> (gamma) e+e-. " \
  "Selection: 2 oppositely charged tracks with p > 2 GeV/c, " \
  "|cos(theta)| < 0.8, EMC energy deposit for saturation-effect handling. " \
  "Cross-checked with di-photon process e+e- -> (gamma) gamma gamma.")

algorithm.note(:bhabha_selection,
  "Bhabha event selection: exactly 2 oppositely charged tracks, " \
  "momentum of each track > 2 GeV/c, |cos(theta)| < 0.8 in EMC barrel. " \
  "EMC energy deposit used to correct saturation effects at high energy. " \
  "These criteria differ from the Lambda_c selection and are applied " \
  "in a separate analysis path.")

algorithm.note(:diphoton_crosscheck,
  "Di-photon cross-check: e+e- -> (gamma) gamma gamma events selected " \
  "with 2 photon candidates in EMC, applying energy and angle cuts. " \
  "Used only for systematic cross-check of luminosity.")

algorithm.with_decay_card(decay_card_lclc).apply(event_selection)
algorithm.execute_on(all_data)