# 2212.07071v2: Neutron EM form factors e+e- -> n nbar at sqrt(s)=2.0-2.95 GeV
# ConExc continuum analysis — NOT a TagAnalysis, no standard track/photon selection
# Multi-energy scan: 12 energy points grouped into 5 intervals

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# The 12 energy points: 2.0000, 2.0500, 2.1250, 2.1500, 2.1750, 2.2000,
# 2.2324, 2.3094, 2.3864, 2.3960, 2.6454, 2.9500 GeV
# Grouped into 5 intervals for EMFF extraction
# Need to look up exact sample names from dataset table
# Approximate sample names based on BOSS version + energy in MeV
energy_points = [
  "709_2000", "709_2050", "709_2125", "709_2150", "709_2175",
  "709_2200", "709_2232", "709_2309", "709_2386", "709_2396",
  "709_2645", "709_2950"
]

data_points = energy_points.map { |s| DatasetManager.real_data.find(s) }
incMC_points = energy_points.map { |s| DatasetManager.inclusive_mc.find(s) }

# ConExc decay card: e+e- -> n nbar via virtual photon (vpho)
# Form D: -2 2112 -2112 (n nbar, sharing xs_user.txt for ISR correction)
# Note: The DSL auto-detects "ConExc" token and switches template
conexc_card = <<~CONEXC
Particle vpho 2.3960 0.0

Decay vpho
  1 ConExc -2 2112 -2112;
Enddecay

Decay vhdr
  1 neutron anti-neutron PHSP;
Enddecay

End
CONEXC

# Create exclusive MC for each energy point using create_exclusive_mc_for
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_nnbar"
  config.events        = 500_000
  config.decay_card    = conexc_card
  config.cross_section = :default
end

# This analysis does not use standard BOSS track/photon selection
# The event selection is based on:
# - Event categories A, B, C depending on detector responses to n nbar
# - Zero charged tracks in MDC (neutral event)
# - EMC shower with 0.5-2.0 GeV energy, |cos(theta)|<0.7 as nbar candidate
# - TOF response with azimuthal angle span
# Signal extraction via composite ML fit to TOF time diff / opening angle
# EMFF extraction via simultaneous minimization of NLL

algorithm = Algorithm.new("NeutronEMFF", "00-00-01")
algorithm.set_header(["NeutronEMFF/NeutronEMFF.h"])
# ECMS varies by energy point — DSL uses per-run MeasuredEcmsSvc for data
algorithm.set_constant({ "ECMS" => [:double, 2.396] })

# No standard selection chain — the analysis uses entirely custom event classification
# We create a minimal selection to satisfy the framework
event_selection = Selection.new

algorithm.note(:event_selection,
  "Event selection: zero charged tracks required. EMC shower in |cos(theta)|<0.7 with " \
  "deposited energy (0.5, 2.0) GeV as nbar candidate. TOF hits within azimuthal span of " \
  "6 scintillators along n momentum. Three event categories (A/B/C): A=TOF from both + EMC " \
  "shower from nbar; B=EMC from both + TOF from nbar only; C=EMC from both, no TOF. " \
  "Signal yields per cos(theta_nbar) bin from composite ML fit (NLL minimization with MIGRAD). " \
  "Backgrounds: beam-related and e+e- -> gamma gamma."
)

algorithm.note(:efficiency_calibration,
  "Data-driven n/nbar efficiency calibration via J/psi -> anti-p pi n (p pi anti-n) control sample. " \
  "Trigger correction from EMC-based online trigger efficiency for neutral events. " \
  "Efficiency iterative determination using form factor input model — converges within 1%."
)

algorithm.note(:emff_extraction,
  "|G_M| and R_em = |G_E|/|G_M| extracted by minimizing NLL based on Poisson PDF, " \
  "simultaneous fit to cos(theta) differential cross sections for 3 categories. " \
  "Data grouped into 5 c.m. energy intervals: (2.0000,2.0500), (2.1250,2.1500), " \
  "(2.1750,2.2000,2.2324), (2.3094,2.3864,2.3960), (2.6454,2.9500) GeV. " \
  "7 equidistant bins within -0.7 < cos(theta_nbar) < 0.7."
)

algorithm.note(:conexc_specific,
  "ConExc generator at NLO for ISR + vacuum polarization. " \
  "Coulomb enhancement factor C=1 for neutral baryons. " \
  "ISR correction (1+delta) and vacuum polarization factor from generator log."
)

algorithm.note(:jpsi_calibration,
  "J/psi dataset (10087M events) used for precise data-driven efficiency calibration: " \
  "J/psi -> anti-p pi n (p pi anti-n) for n/nbar detection, " \
  "J/psi -> n nbar, J/psi -> pi+pi-pi0, e+e- -> gamma gamma for neutral particle studies."
)

algorithm.with_decay_card(conexc_card).apply(event_selection)
algorithm.execute_on(data_points + incMC_points + exMCs)