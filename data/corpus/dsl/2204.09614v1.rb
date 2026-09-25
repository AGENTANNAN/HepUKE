### Dataset description ###
# Paper: Amplitude analysis and branching fraction measurement of D_s+ -> K_S0 K+ pi0
# arXiv: 2204.09614v1
# Uses single-tag (ST) and double-tag (DT) with D_s tags at sqrt(s)=4.178-4.226 GeV
# 6.32 fb^{-1} total integrated luminosity

data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4226 = DatasetManager.real_data.find("703_4226")
all_data = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4226]

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4226")
all_incMC = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4226]

# e+e- -> D_s*+ D_s- -> gamma D_s+ D_s-
# Tag D_s- via hadronic modes; signal D_s+ -> K_S0 K+ pi0
# The transition photon gamma from D_s*+ decay is part of the tag reconstruction
decay_card_ds = <<~DECAYCARD
    Decay D_s*+
    1.0000 D_s+ gamma VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.0000 K_S0 K+ pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC generated at a representative energy
exMC_ds = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_dsp_kskpipi0"
  config.related_dataset = data_4226    # primary energy point
  config.events = 2000000
  config.decay_card = decay_card_ds
  config.cross_section = :default
end

### Event selection (BOSS) — TagAnalysis ###
alg_name = "DsTagKsKPi0"
my_algorithm = TagAnalysis.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 4.226] })

# Tag side: D_s- reconstructed via 7 hadronic modes (authoritative list)
# Paper uses 8 modes; the 8th (K_S0 K+ pi- pi-) is not available in the frozen release
my_algorithm.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKKPiPi0, :DstoKsKPiPi,
          :DstoPiPiPi, :DstoPiEtaPrime, :DstoKPiPi
  t.charm -1
end

# Signal side: D_s+ -> K_S0 K+ pi0
# K_S0 -> pi+ pi-: two charged pions tracked at signal side
# pi0 -> gamma gamma: two photons
# K+: one charged kaon
my_algorithm.signal_side do |s|
  s.photons 2
  s.charged(kp: 1, pip: 1, pim: 1)
  s.require_charge 1          # D_s+ charge +1 = K+ + pi+ + pi- = +1
end

my_algorithm.note(:ks_reconstruction, "K_S0 reconstructed via pi+pi- vertex fit; pions must satisfy |Vz|<20 cm, |cos_theta|<0.93, vertex fit chi2<100; K_S0 decay length > 2 sigma; invariant mass in (0.492,0.503) GeV/c2. NOTE: TagAnalysis DSL v1 cannot express K_S0 vertex fit on signal side — this is handled by DTagTool for tag-side K_S0 modes, but signal-side K_S0 reconstruction requires custom code")
            .note(:pi0_eta_reconstruction, "pi0/eta -> gamma gamma: E_gamma > 25 MeV (barrel, |cos_theta|<0.80), > 50 MeV (endcap, 0.86<|cos_theta|<0.92); timing < 700 ns. 1C kinematic fit constraining M(gamma gamma) to known pi0/eta mass, chi2 < 10")
            .note(:eta_prime_reconstruction, "eta' -> pi+ pi- eta with eta -> gamma gamma; M(pi+pi-eta) within 10 MeV of known eta' mass")
            .note(:transition_photon, "Transition photon gamma from D_s*+ -> gamma D_s+ is reconstructed; photon energy requirement applied")
            .note(:tag_selection, "ST candidate: M_tag in [1.930, 1.990] GeV/c2 for D_s-; best ST candidate per tag mode chosen by recoil mass closest to known D_s*+ mass")
            .note(:dt_selection, "DT candidate: best chosen by average mass (M_tag + M_sig)/2 closest to known D_s mass per tag mode")
            .note(:kinematic_fit, "8C kinematic fit: 4-momentum conservation + mass constraints on K_S0, pi0, D_s-, D_s*+; chi2 minimum candidate chosen; then 9th constraint on D_s+ mass added")
            .note(:dt_signal_region, "M_sig in (1.930, 1.990) GeV/c2 for accepted DT candidates")
            .note(:combinatorial_suppression, "Intermediate resonance cuts: M(K+K-)<1.05 GeV (phi), |M(K+pi-)-M(K*(892))|<70 MeV (K*892), |M(pi-pi0)-M(rho)|<150 MeV (rho)")
            .note(:ds_kpi_purity, "For D_s+ -> K_S0 K+ pi0 signal: D_s+ momentum > 0.1 GeV/c to remove soft pions from D*+ decays")
            .note(:missing_tag_mode, "The 8th tag mode K_S0 K+ pi- pi- (D_s- decay) is not in the frozen BOSS 7.0.6/7.1.2 DTagAlg enum and cannot be expressed in TagAnalysis DSL v1")

my_algorithm.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

my_algorithm.apply
root_files = my_algorithm.execute_on(all_data + all_incMC + [exMC_ds])