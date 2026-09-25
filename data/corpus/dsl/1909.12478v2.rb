### Dataset preparation ###
# Seven energy points for e+e- -> D+ D- pi+ pi- cross-section measurement
# Energy points (MeV): 4358.3, 4387.4, 4415.6, 4467.1, 4527.1, 4574.5, 4599.5
# BOSS version 703 for psi(4260) energy region scans

datasets = [
  DatasetManager.real_data.find("703_4358"),
  DatasetManager.real_data.find("703_4387"),
  DatasetManager.real_data.find("703_4416"),
  DatasetManager.real_data.find("703_4467"),
  DatasetManager.real_data.find("703_4527"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
]

incMC_datasets = [
  DatasetManager.inclusive_mc.find("703_4358"),
  DatasetManager.inclusive_mc.find("703_4387"),
  DatasetManager.inclusive_mc.find("703_4416"),
  DatasetManager.inclusive_mc.find("703_4467"),
  DatasetManager.inclusive_mc.find("703_4527"),
  DatasetManager.inclusive_mc.find("703_4575"),
  DatasetManager.inclusive_mc.find("703_4600"),
]

# Decay card: e+e- -> D+ D- pi+ pi- with intermediate D1(2420)+ -> D+ pi+ pi-
# Top mother: psi(4260) (KKMC convention for continuum region)
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 D+ D- pi+ pi- PHSP;
    Enddecay
    Decay D+
    1.000 K- pi+ pi+ PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC for multiple energy points (create_exclusive_mc_for)
exMCs = DatasetManager.create_exclusive_mc_for(datasets) do |config|
  config.sample_name = "ee_to_Dp_Dm_pi_pi"
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection for e+e- -> D+ D- pi+ pi- ###

alg = Algorithm.new("EEToDpDmPiPi")
alg.set_header(["EEToDpDmPiPiAlg/EEToDpDmPiPi.h"])
   .set_constant({"ECMS" => [:double, 4.416]})   # Nominal value; per-run ECMS from MeasuredEcmsSvc
   .note(:multi_energy, "Analysis covers 7 energy points: 4358.3, 4387.4, 4415.6, 4467.1, 4527.1, 4574.5, 4599.5 MeV with per-run ECMS")
   .note(:recoil_mass_analysis, "Recoil mass RM(D+ pi+ pi-) used to identify D- signal; RM(D+) distribution fitted for D1(2420)+ and psi(3770)pi+pi- signals; sideband subtraction for combinatorial background. Recoil-mass technique not directly expressible in DSL")
   .note(:dplus_vertex_fit, "Vertex fit for D+ -> K- pi+ pi+ candidates (chi2_VF < 100) followed by kinematic fit constraining D+ mass to nominal (chi2_KF < 20)")
   .note(:cross_section, "Born cross sections measured for D1(2420)+D-+c.c. and psi(3770)pi+pi- processes; radiative correction factors (1+delta_rad) and vacuum polarization (1/|1-Pi|^2) applied in ROOT analysis")

sel = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=3"
    nChrn ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkm ">=1"
  }
  .remove([:km <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  # Kinematic fit for D+ -> K- pi+ pi+ mass constraint (vertex fit implied by vertex_fit)
  .kinematic_fit([:km, :pip, :pip]) {
    nominal
    invariant_mass_of(:km, :pip, :pip).constrain_to_nominal_mass_of(:Dp)
    constrain_four_momentum
    chi2_cut 20
  }

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on(datasets + incMC_datasets + exMCs)