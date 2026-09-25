# BESIII measurement of nbar p -> K+K-pi+ (pi0) cross sections using antineutrons
# from J/psi -> p pi- nbar decays. Data: J/psi sample at sqrt(s)=3.097 GeV.

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---------------- Decay cards ----------------
# Signal: J/psi -> p pi- nbar (antineutron tagging), followed by nbar p scattering
# on hydrogen in oil layer producing K+ K- pi+ (pi0).
# The scattering is generated in phase space.

dc_KKpi = <<~DC
  Decay J/psi
  1.0000 p+ pi- anti-n0 PHSP;
  Enddecay
  End
DC

dc_KKpipi0 = dc_KKpi

exmc_tag = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Jpsi_ppinbar_tag"
  c.related_dataset = jpsi_data
  c.events          = 2_000_000
  c.decay_card      = dc_KKpi
  c.cross_section   = :default
end

# ============ nbar p -> K+ K- pi+ ============
alg1 = Algorithm.new("NbarPtoKKPi")
alg1.set_header(["NbarPtoKKPiAlg/NbarPtoKKPi.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })

sel1 = Selection.new
sel1.select_track {
       cos_theta 0.93
       Vz 10.0
       Vr 1.0
       nChrp ">=2"
       nChrn ">=2"
     }
    .pid(method: :probability) {
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       identify :kaon,   against: [:proton, :pion]
       identify :pion,   against: [:proton, :kaon]
       nprp ">=1"
       npim ">=1"
       nkp ">=1"
       nkm ">=1"
       npip ">=1"
     }
    .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
     }
    .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:p_pi_tag).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
    # 1C kinematic fit constraining the recoil against p pi- to the antineutron mass
    .kinematic_fit([:prp, :pim, :kp, :km, :pip]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
     }

alg1.note(:antineutron_tagging,
          "tag antineutron via J/psi -> p pi-; vertex fit on p pi- with chi2 < 5, then 1C kinematic " \
          "fit constraining RM(p pi-) to nominal nbar mass; require kinematic-fit chi2 < 11")
    .note(:secondary_vertex_kkpi,
          "vertex fit of K+K-pi+ (chi2 < 6); vertex position taken as nbar-p scattering point (x1,y1,z1)")
    .note(:scattering_position,
          "require |Delta r| < 0.3 cm and |Delta z| < 2.0 cm between fitted vertex and expected " \
          "antineutron impact point on the beam pipe")
    .note(:proton_momentum_cut,
          "P(p) = |P(K+K-pi+) - P(nbar)| < 0.04 GeV/c to select scattering on 1H in the oil layer")
alg1.with_decay_card(dc_KKpi).apply(sel1)
alg1.execute_on([jpsi_data, jpsi_incMC, exmc_tag])

# ============ nbar p -> K+ K- pi+ pi0 ============
alg2 = Algorithm.new("NbarPtoKKPiPi0")
alg2.set_header(["NbarPtoKKPiPi0Alg/NbarPtoKKPiPi0.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })

sel2 = Selection.new
sel2.select_track {
       cos_theta 0.93
       Vz 10.0
       Vr 1.0
       nChrp ">=2"
       nChrn ">=2"
     }
    .pid(method: :probability) {
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       identify :kaon,   against: [:proton, :pion]
       identify :pion,   against: [:proton, :kaon]
       nprp ">=1"
       npim ">=1"
       nkp ">=1"
       nkm ">=1"
       npip ">=1"
     }
    .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam ">=2"
     }
    .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 50
       npi0 ">=1"
     }
    .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:p_pi_tag).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
    .kinematic_fit([:prp, :pim, :kp, :km, :pip, :pi0]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
     }

alg2.note(:antineutron_tagging,
          "tag antineutron via J/psi -> p pi- nbar; vertex fit on p pi- with chi2 < 5, then 1C fit " \
          "constraining RM(p pi-) to nbar mass; require kinematic-fit chi2 < 11")
    .note(:secondary_vertex_kkpipi0,
          "vertex fit of K+K-pi+ (chi2 < 6); vertex position taken as nbar-p scattering point")
    .note(:scattering_position,
          "require |Delta r| < 0.3 cm and |Delta z| < 2.0 cm between fitted vertex and expected " \
          "antineutron impact point on the beam pipe")
    .note(:proton_momentum_cut,
          "P(p) = |P(K+K-pi+ pi0) - P(nbar)| < 0.04 GeV/c to select scattering on 1H in the oil layer")
alg2.with_decay_card(dc_KKpipi0).apply(sel2)
alg2.execute_on([jpsi_data, jpsi_incMC, exmc_tag])
