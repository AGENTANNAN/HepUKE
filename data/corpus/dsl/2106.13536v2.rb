DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
all_data = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
all_incMC = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0 D_s*+ D_s- PHSP;
  Enddecay
  Decay D_s*+
  1.0 gamma D_s+ VSP_PWAVE 1.0 2.0 0.0;
  Enddecay
  Decay D_s+
  1.0 pi+ pi+ pi- eta PHSP;
  Enddecay
  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_Ds_3pipi_eta"
  config.events        = 500_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

algorithm = TagAnalysis.new("DsTag3PiPiEta", "00-00-01")
algorithm
  .set_header(["DsTag3PiPiEta/DsTag3PiPiEta.h"])
  .set_constant({ "ECMS" => [:double, 4.200] })
  .note(:analysis_overview,
    "First observation of Ds+ -> pi+ pi+ pi- eta and amplitude analysis " \
    "using 6.32 fb^-1 at ECM = 4.178-4.226 GeV. " \
    "Double-tag technique: tag Ds- from 8 hadronic modes, " \
    "signal Ds+ -> pi+ pi+ pi- eta reconstructed from remaining tracks. " \
    "Measured BF = (3.12 +/- 0.13 +/- 0.09)%. " \
    "Observed W-annihilation Ds+ -> a0(980)+ rho(770)0 with BF (0.21 +/- 0.08 +/- 0.05)%.")
  .note(:tag_modes,
    "Eight hadronic tag modes for Ds-: " \
    "DstoKsK, DstoKKPi, DstoKsKPi0, DstoKKPiPi0, " \
    "DstoKsKminusPiPi, DstoKsKplusPiPi, DstoPiPiPiEta, DstoKPiPi. " \
    "For Ds- -> K- pi+ pi- mode, M(pi+pi-) in [0.487,0.511] GeV/c^2 excluded " \
    "to avoid overlap with Ds- -> K_S0 K-.")
  .note(:kinematic_fit,
    "7C kinematic fit: 4-momentum conservation + eta mass constraint + " \
    "tag Ds- mass constraint + Ds*+ mass constraint. " \
    "D_s*+ mass constraint involving transition photon is applied in ROOT; " \
    "DSL fit constrains 4-momentum, eta, and tag Ds mass. " \
    "Multiple candidates (approx 15%): select minimum chi2. " \
    "For BF measurement: transition photon from Ds* not reconstructed, " \
    "BDT requirement dropped to improve statistical precision.")
  .note(:background_suppression,
    "K_S0 veto: secondary vertex fit on pi+pi- pair, " \
    "reject if M_pipi in [0.487,0.511] GeV/c^2 and L/sigma_L > 2. " \
    "eta' veto: reject M(pi+pi-eta) < 1 GeV/c^2. " \
    "pi0 cross-feed veto: reject if M(gamma_eta gamma_pi0) or " \
    "M(gamma_eta gamma_other) in [0.115,0.150] GeV/c^2. " \
    "BDT with 4 input variables to ensure purity > 85%. " \
    "ST yield = 479,093 +/- 1,952; DT signal yield = 2,139 +/- 78.")
  .note(:amplitude_analysis,
    "Unbinned ML amplitude analysis with 11 intermediate amplitudes. " \
    "Dominant: Ds+ -> a1(1260)+ eta, a1(1260)+ -> rho(770)0 pi+ (FF=55.4%). " \
    "W-annihilation: Ds+ -> a0(980)+ rho(770)0 (FF=6.7%). " \
    "Goodness of fit: chi2/NDOF = 153.2/133 = 1.15, p-value = 11.1%.")
  .with_decay_card(decay_card)

algorithm.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0, :DstoKKPiPi0,
          :DstoKsKminusPiPi, :DstoKsKplusPiPi, :DstoPiPiPiEta, :DstoKPiPi
  t.charm -1
end

algorithm.signal_side do |s|
  s.photons 2
  s.charged(pip: 2, pim: 1)
  s.min_photon_angle 10.0
end

algorithm.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D_s+")
  f.chi2_cut 200
end

algorithm.apply
algorithm.execute_on(all_data + all_incMC + [sig_mc])