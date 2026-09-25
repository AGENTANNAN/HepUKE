# 2210.09561v2: Inclusive Λc- → nbar + X
# ST method: tag Λc+ → pK-π+, signal anti-neutron via EMC shower
# 7 energy points: 4.600-4.699 GeV
# Anti-neutron selection via shower properties: E_nbar>0.48GeV, H_nbar>20, S_nbar>18cm²

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Data at 7 c.m. energies (Table I)
data_4600 = DatasetManager.real_data.find("703_4600")
data_4612 = DatasetManager.real_data.find("706_4610")
data_4628 = DatasetManager.real_data.find("706_4620")
data_4641 = DatasetManager.real_data.find("706_4640")
data_4661 = DatasetManager.real_data.find("706_4660")
data_4682 = DatasetManager.real_data.find("706_4680")
data_4699 = DatasetManager.real_data.find("706_4700")

all_data = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699]

alg = TagAnalysis.new("LcInclusiveNbar")
alg.set_header(["LcInclusiveNbarAlg/LcInclusiveNbar.h"])
   .set_constant({ "ECMS" => [:double, 4.640] })

# ST: Λc+ → pK-π+ (one tag mode)
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP
  t.charm 1
  t.window :deltaE, min: -0.034, max: 0.020
end

# Signal side: no proton tracks allowed; anti-neutron detected via EMC shower
alg.signal_side do |s|
  s.photons 0
  s.charged(prp: 0, prm: 0, at_least: true)
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply

# ======================================================================
# EXTENSIVE NOTES: The bulk of the anti-neutron selection is inexpressible
# in the current TagAnalysis DSL because it involves shower-shape variables
# and data-driven efficiency corrections that BOSS does not handle natively.
# ======================================================================
alg.note(:nbar_shower_selection,
  "Anti-neutron identification via EMC shower: " \
  "most energetic shower required to have " \
  "E_nbar > 0.48 GeV, H_nbar > 20 (hit crystals), " \
  "S_nbar > 18 cm² (second moment of shower shape). " \
  "Angle between charged track and shower > 20°. " \
  "These cuts are applied in the generated C++ code by hand.")

alg.note(:nbar_data_driven_efficiency,
  "Detection efficiency of anti-neutron selections corrected via a " \
  "data-driven method using a control sample of J/ψ → p nbar π- at " \
  "√s = 3.097 GeV. Efficiency applied as weight in bins of " \
  "anti-neutron momentum and cosθ.")

alg.note(:no_extra_proton,
  "Signal side requires no tracks identified as proton (PID highest " \
  "probability) within Vz < 20 cm to suppress decays with a proton " \
  "in the final state.")

alg.note(:st_selection,
  "ST Λc+ candidates selected via MBC ∈ (2.275, 2.31) GeV/c² and " \
  "ΔE ∈ (−34, 20) MeV (asymmetric, 3σ). If multiple candidates, " \
  "keep the one with minimum |ΔE|.")

alg.note(:background_subtraction,
  "qqbar background estimated from MBC sideband (2.20, 2.26) GeV/c². " \
  "Λc+Λc- background from inclusive MC of Λc- → pbar + X. " \
  "Signal yield from MBC fit after applying nbar selections.")

alg.execute_on(all_data)