# Observation of near-threshold enhancement in Lambda anti-Lambda
# from e+e- -> phi Lambda anti-Lambda at sqrt(s) 3.51-4.60 GeV
# Paper: 2104.08754v2
# Multi-energy scan at 28 energy points
# Partial reconstruction: one Lambda/anti-Lambda missing,
# kinematic fit constrains missing mass to Lambda nominal mass

algo = Algorithm.new('ee_to_phi_Lambda_antiLambda', version: '00-00-01')
algo.set_header(['VertexFit/VertexDbSvc.h'])
algo.set_constant(ECMS: 4.178)

algo.with_decay_card <<~DECAY
  Decay vpho
  1.000 phi Lambda anti-Lambda PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- PHSP;
  Enddecay
  Decay Lambda
  1.000 p+ pi- PHSP;
  Enddecay
  Decay anti-Lambda
  1.000 anti-p- pi+ PHSP;
  Enddecay
DECAY

event_selection = Selection.new('PhiLambdaAntiLambda') do |sel|
  # 4 charged tracks: p, pi-, K+, K- (one Lambda missing)
  sel.select_track do |t|
    t.nTot = 4
    t.nChrp = 2
    t.nChrn = 2
    t.cos_theta_max = 0.93
    t.vr_max = 10.0
    t.vz_max = 20.0
  end

  # Probability-based PID for all hadrons
  sel.pid(method: :probability) do |pid|
    pid.identify(:prp, :prm, against: [:pip, :pim, :kp, :km])
    pid.identify(:kp, :km, against: [:pip, :pim, :prp, :prm])
    pid.prob_cut(0.001)
  end

  # Reconstruct Lambda -> p+ pi- with secondary vertex fit
  sel.secondary_vertex_fit(:Lambda, daughters: [:prp, :pim]) do |vtx|
    vtx.chi2_cut(100)
  end

  # Partial reconstruction: missing anti-Lambda constrained to Lambda mass
  # reconstruct e+e- -> phi + Lambda(seen) + anti-Lambda(missing)
  sel.partial_rec do |pr|
    pr.miss_particle(:Lambda_bar)
    pr.seen_particles([:Lambda, :kp, :km])
    pr.chi2_cut(30)
  end
end

algo.apply(event_selection)

# Exclusive MC (at primary energy point for efficiency study)
DatasetManager.create_exclusive_mc do |mc|
  mc.decay_card <<~DECAY
    Decay vpho
    1.000 phi Lambda anti-Lambda PHSP;
    Enddecay
    Decay phi
    1.000 K+ K- PHSP;
    Enddecay
    Decay Lambda
    1.000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda
    1.000 anti-p- pi+ PHSP;
    Enddecay
  DECAY
  mc.related_dataset DatasetManager.load_real_data.find('703_4180')
  mc.n_events 500_000
end

# Execute on 28 energy points from 3.51 to 4.60 GeV
algo.execute_on([
  # BOSS 703 data
  DatasetManager.load_real_data.find('703_3510'),  # 3.5106 GeV (chi_c1_scan_4)
  DatasetManager.load_real_data.find('712_3773'),  # 3.773 GeV (psi(3770))
  DatasetManager.load_real_data.find('703_3872'),  # 3.867 GeV
  DatasetManager.load_real_data.find('703_4009'),  # 4.0076 GeV
  # BOSS 705 data
  DatasetManager.load_real_data.find('705_4130'),  # 4.1285 GeV
  DatasetManager.load_real_data.find('705_4160'),  # 4.1574 GeV
  # BOSS 703 data
  DatasetManager.load_real_data.find('703_4180'),  # 4.178 GeV
  DatasetManager.load_real_data.find('703_4190'),  # 4.1888 GeV
  DatasetManager.load_real_data.find('703_4200'),  # 4.1989 GeV
  DatasetManager.load_real_data.find('703_4210'),  # 4.2092 GeV
  DatasetManager.load_real_data.find('703_4220'),  # 4.2187 GeV
  DatasetManager.load_real_data.find('703_4230'),  # 4.2263 GeV
  DatasetManager.load_real_data.find('703_4237'),  # 4.2357 GeV
  DatasetManager.load_real_data.find('703_4246'),  # 4.2438 GeV
  DatasetManager.load_real_data.find('703_4260'),  # 4.2580 GeV
  DatasetManager.load_real_data.find('703_4270'),  # 4.2668 GeV
  DatasetManager.load_real_data.find('703_4280'),  # 4.2777 GeV
  # BOSS 705 data
  DatasetManager.load_real_data.find('705_4290'),  # 4.2879 GeV
  DatasetManager.load_real_data.find('705_4315'),  # 4.3120 GeV
  DatasetManager.load_real_data.find('705_4340'),  # 4.3374 GeV
  # BOSS 703 data
  DatasetManager.load_real_data.find('703_4360'),  # 4.3583 GeV
  # BOSS 705 data
  DatasetManager.load_real_data.find('705_4380'),  # 4.3774 GeV
  DatasetManager.load_real_data.find('705_4400'),  # 4.3965 GeV
  # BOSS 703 data
  DatasetManager.load_real_data.find('703_4420'),  # 4.4156 GeV
  # BOSS 705 data
  DatasetManager.load_real_data.find('705_4440'),  # 4.4362 GeV
  # BOSS 703 data
  DatasetManager.load_real_data.find('703_4470'),  # 4.4671 GeV
  DatasetManager.load_real_data.find('703_4530'),  # 4.5271 GeV
  DatasetManager.load_real_data.find('703_4600'),  # 4.5995 GeV
  # Inclusive MC
  DatasetManager.load_inclusive_mc.find('703_4180')
])

algo.note(:partial_reconstruction,
  'One Lambda/anti-Lambda is assumed missing. Kinematic fit constrains ' \
  'the missing mass to the nominal Lambda mass. ' \
  'Events with multiple candidates: choose combination with smallest ' \
  'combined chi2 from vertex and kinematic fits.')
algo.note(:vertex_fit,
  'Lambda secondary vertex fit: p pi- from common vertex. ' \
  'Second vertex fit for K+ K- Lambda from IP. ' \
  'Flight distance significance > 2 sigma between IP and Lambda decay vertex.')
algo.note(:mass_windows,
  'Lambda mass window: [1.112, 1.120] GeV/c^2 (about +/- 3 sigma). ' \
  'phi mass window: [1.01, 1.03] GeV/c^2. ' \
  'Combined vertex + kinematic fit chi2 < 30.')
algo.note(:multi_energy,
  'Combined analysis at 28 energy points from 3.51 to 4.60 GeV. ' \
  'Total integrated luminosity: 19.5 fb^-1. ' \
  'Cross section measurement and lineshape analysis performed in ROOT stage.')
algo.note(:signal_extraction,
  'Signal yields extracted by fitting K+K- invariant mass spectrum ' \
  'with MC-simulated signal shape + inverted ARGUS background. ' \
  'Lambda anti-Lambda near-threshold enhancement studied via ' \
  'Breit-Wigner and reversed exponential fits in ROOT analysis.')