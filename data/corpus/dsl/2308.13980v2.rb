# Paper 2308.13980v2: Search for χc1(3872)→π+π-η at √s 4.13-4.34 GeV
# Two η decay modes: η→γγ (Mode I) and η→π+π-π0 (Mode II)
# Multi-energy scan, 16 energy points, 11.5 fb-1

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Energy points from the paper Table I (BOSS 703/705)
data_points = [
  DatasetManager.real_data.find("703_4130"),  # 4128.8 MeV
  DatasetManager.real_data.find("703_4160"),  # 4157.8 MeV
  DatasetManager.real_data.find("703_4180"),  # 4178.0 MeV
  DatasetManager.real_data.find("703_4190"),  # 4189.1 MeV (526.7/pb)
  DatasetManager.real_data.find("703_4200"),  # 4199.2 MeV
  DatasetManager.real_data.find("703_4210"),  # 4209.4 MeV (517.1/pb)
  DatasetManager.real_data.find("703_4220"),  # 4218.9 MeV
  DatasetManager.real_data.find("703_4230"),  # 4226.3 MeV
  DatasetManager.real_data.find("703_4237"),  # 4235.8 MeV
  DatasetManager.real_data.find("703_4246"),  # 4244.0 MeV
  DatasetManager.real_data.find("703_4260"),  # 4258.0 MeV
  DatasetManager.real_data.find("703_4270"),  # 4266.8 MeV
  DatasetManager.real_data.find("703_4280"),  # 4277.8 MeV
  DatasetManager.real_data.find("705_4290"),  # 4288.4 MeV
  DatasetManager.real_data.find("705_4315"),  # 4312.7 MeV
  DatasetManager.real_data.find("705_4340"),  # 4337.9 MeV
]

inc_mc_points = data_points.map { |d|
  begin
    DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
  rescue
    nil
  end
}.compact

# ============================================================
# Mode I: η → γγ
# Final state: γ_rad + π+ + π- + γ + γ
# 5C kinematic fit (4C + η mass constraint)
# ============================================================

decay_card_mode1 = <<~DECAYCARD
  Decay vpho
  1.0000 gamma chi_c1(3872)  VSP_PWAVE 1 0 0 1;
  Enddecay
  Decay chi_c1(3872)
  1.0000 pi+ pi- eta  PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

algo_mode1 = Algorithm.new("x3872_pipi_eta_gg", "00-00-01")
  .set_header(["x3872_pipi_eta_gg/x3872_pipi_eta_gg.h"])
  .with_decay_card(decay_card_mode1)

exMC_mode1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_x3872_pipi_eta_gg"
  config.events        = 100_000
  config.decay_card    = decay_card_mode1
  config.cross_section = :default
end

selection_mode1 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .assign({chrgp: :pip, chrgn: :pim})
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  end
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_mode1
  .note(:pi0_veto, "|M(γ_rad γ_1(2)) - m(π0)| > 20 MeV/c^2 to suppress π0→γγ mis-combination backgrounds")
  .note(:dimuon_veto, "cosθ_ππ > -0.96 to reject radiative dimuon background in η→γγ mode")
  .apply(selection_mode1)

# ============================================================
# Mode II: η → π+π-π0, π0 → γγ
# Final state: γ_rad + π+ + π- + π+ + π- + γ + γ
# 5C kinematic fit (4C + η mass constraint)
# ============================================================

decay_card_mode2 = <<~DECAYCARD
  Decay vpho
  1.0000 gamma chi_c1(3872)  VSP_PWAVE 1 0 0 1;
  Enddecay
  Decay chi_c1(3872)
  1.0000 pi+ pi- eta  PHSP;
  Enddecay
  Decay eta
  1.0000 pi+ pi- pi0  DALITZ;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

algo_mode2 = Algorithm.new("x3872_pipi_eta_3pi", "00-00-01")
  .set_header(["x3872_pipi_eta_3pi/x3872_pipi_eta_3pi.h"])
  .with_decay_card(decay_card_mode2)

exMC_mode2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_x3872_pipi_eta_3pi"
  config.events        = 100_000
  config.decay_card    = decay_card_mode2
  config.cross_section = :default
end

selection_mode2 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
  end
  .pid(method: :probability) do
    identify :pion, against: [:kaon]
    npip ">=1"  # at least one pi+ identified as pion vs kaon
    npim ">=1"  # at least one pi- identified as pion vs kaon
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :gamma]) do
    invariant_mass_of(:pip, :pim, :pip, :pim, :pi0).constrain_to_nominal_mass_of(:eta)
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_mode2
  .note(:radiative_photon_veto, "veto events with M(γ_rad γ_1) in [125.6,150.0] and M(γ_rad γ_2) in [115.7,160.0] MeV/c^2 to suppress e+e-→π+π-π+π-π0 background")
  .apply(selection_mode2)

# Execute on datasets
all_exMC = exMC_mode1 + exMC_mode2

algo_mode1.execute_on(data_points + inc_mc_points + exMC_mode1)
algo_mode2.execute_on(data_points + inc_mc_points + exMC_mode2)