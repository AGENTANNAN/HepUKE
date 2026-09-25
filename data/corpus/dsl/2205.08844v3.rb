# 2205.08844v3: Amplitude analysis and branching fraction measurement
# of Ds+ -> K+ pi+ pi-
#
# Tag-based (Ds double-tag). Uses 6.32 fb-1 at sqrt(s)=4.178-4.226 GeV.
# ST: Ds- in 10 hadronic tag modes. DT: tag + signal Ds+ -> K+ pi+ pi-.
# Ds* -> gamma Ds transition photon included. 6C/7C kinematic fit.
# KS0 veto on pi+pi- mass.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# --- Datasets (classified into 3 sample groups per paper) ---
data_4178 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4226 = DatasetManager.real_data.find("703_4230")

incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190scan") rescue DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210scan") rescue DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220scan") rescue DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

all_data    = [data_4178, data_4190, data_4200, data_4210, data_4220, data_4226]
all_incMC   = [incMC_4178, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4226]

# --- Decay card for e+e- -> Ds*+ Ds- -> gamma Ds+ Ds- ---
decay_card_ds = <<~DECAYCARD
  Decay psi(4260)
  1.0 D_s*+ D_s- PHSP;
  Enddecay
  Decay D_s*+
  1.0 gamma D_s+ VSP_PWAVE;
  Enddecay
  Decay D_s*-
  1.0 gamma D_s- VSP_PWAVE;
  Enddecay
  End
DECAYCARD

# --- TagAnalysis: Ds DT ---
alg = TagAnalysis.new("DsDstTagKPiPi")
alg.set_header(["DsDstTagKPiPiAlg/DsDstTagKPiPi.h"])
    .set_constant({ "ECMS" => [:double, 4.200] })
    .with_decay_card(decay_card_ds)

# Tag side 1: Ds- in 10 hadronic tag modes
alg.tag_side(:Ds) do |t|
  t.modes :DstoKPi0, :DstoKKPi, :DstoKPi0Pi0,
          :DstoKKPiPi0, :DstoKKPiPiPi, :DstoKsK,
          :DstoPiPiPi, :DstoPiEta, :DstoPiEtaPrime,
          :DstoKSEta
  t.charm -1
end

# Tag side 2: Ds+ in the signal mode Ds+ -> K+ pi+ pi-
alg.tag_side(:Ds) do |t|
  t.modes :DstoKPiPi
  t.charm 1
  t.rank_by :inv
end

# Signal side: transition photon from Ds* -> gamma Ds
alg.signal_side do |s|
  s.photons 1
  s.min_photon_angle 10.0
end

# Fit: 6C kinematic fit with Ds* mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag1, :gamma).constrain_to_nominal_mass_of(:"D_s*+")
  f.invariant_mass_of(:tag2, :gamma).constrain_to_nominal_mass_of(:"D_s*+")
  f.chi2_cut 200
end

alg.note(:ks0_veto,
  "K_S0 veto: pi+ pi- invariant mass outside [0.4676, 0.5276] GeV/c2 " \
  "to reject Ds+ -> KS0 K+ background. Applied in ROOT stage after the " \
  "7C (6C + signal Ds mass constraint) kinematic fit used for amplitude analysis.")

alg.note(:mrec_windows,
  "M_rec windows per energy point used to suppress non-Ds*Ds background: " \
  "[2.050,2.180] @ 4.178, [2.048,2.190] @ 4.189, [2.046,2.200] @ 4.199, " \
  "[2.044,2.210] @ 4.209, [2.042,2.220] @ 4.219, [2.040,2.220] @ 4.226 GeV/c2. " \
  "Applied in ROOT stage.")

alg.note(:amplitude_analysis,
  "Amplitude analysis via unbinned maximum likelihood fit in ROOT. " \
  "Three sample groups: 4.178, 4.189-4.219, and 4.226 GeV. " \
  "Signal yields: 772, 444, 140 with purities ~93-96%. " \
  "BF = (6.11 +/- 0.18_stat +/- 0.11_syst) x 10^-3.")

alg.apply
alg.execute_on(all_data + all_incMC)