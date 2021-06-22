from pyrosetta import *
from pyrosetta.rosetta import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--pdb', '-p', type=str, required=True, help='Starting PDB model')
parser.add_argument('--resfile', '-r', type=str, required=True, help='resfile path')
parser.add_argument('--symdef', type=str, required=True, help='symdef file path')
parser.add_argument('--fc', type=str, default=False, help='ignore first 207 residues? True/False')
parser.add_argument('--disallow_resis', '-d', type=str, help='disallow any residues, e.g. CPGWYH are common. Default off.')
parser.add_argument('--ld', '-l', type=str, default='True', help='Turn layer design on? Accepts on or off')
args = parser.parse_args()

pyrosetta.init("-beta -unmute all -ex1 -ex2aro -overwrite -native " + args.pdb + " -s " + args.pdb)

pose = pose_from_file(args.pdb)
output_name = args.pdb.replace('.pdb', '_0001.pdb.gz')
#resfile = './resfile.txt'
fc_on = args.fc

scorefxn = core.scoring.ScoreFunctionFactory.create_score_function('beta')
sym_scorefxn = core.scoring.symmetry.SymmetricScoreFunction()
sym_scorefxn.assign(scorefxn)

def set_up_tf(fc_on=False):
	
	limit_aro = protocols.task_operations.LimitAromaChi2Operation()
	limit_aro.chi2max(110)
	limit_aro.chi2min(70)

	lock_PG = core.select.residue_selector.ResidueNameSelector()
	lock_PG.set_residue_names("PRO,GLY")

	#disallowed = core.pack.task.operation.DisallowIfNonnative()
	#disallowed.disallow_aas('CPGW')

	tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
	tf.push_back(core.pack.task.operation.InitializeFromCommandline())
	tf.push_back(core.pack.task.operation.IncludeCurrent())
	tf.push_back(limit_aro)
	if args.disallow_resis:
		#disallowed = ''
		#for resi in args.disallow_resis:
		#	disallowed += resi
		disallow_if_nonnative = core.pack.task.operation.DisallowIfNonnative()
		#disallow_if_nonnative.disallow_aas(disallowed)
		disallow_if_nonnative.disallow_aas(args.disallow_resis)
		tf.push_back(disallow_if_nonnative)
	tf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(args.resfile))
	if fc_on == True:
		fc_range = core.select.residue_selector.ResidueIndexSelector('1-207')
		#no_fc = core.select.residue_selector.NotResidueSelector(fc_range)
		tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(\
		pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), fc_range))
	tf.push_back(core.pack.task.operation.OperateOnResidueSubset(\
		core.pack.task.operation.PreventRepackingRLT(), lock_PG))
	return(tf)

def add_layer_design_to_tf():
	task_f = set_up_tf()

	core_resis = core.select.residue_selector.LayerSelector()
	core_resis.set_layers(True, False, False) #sets core
	core_resis.use_sc_neighbors()
	core_resis.set_use_sc_neighbors(True)
	core_resis.set_cutoffs(4.9, 2.7)
	core_resis.set_dist_exponent(0.7)

	boundary_resis = core.select.residue_selector.LayerSelector()
	boundary_resis.set_layers(False, True, False) #sets boundary
	boundary_resis.use_sc_neighbors()
	boundary_resis.set_use_sc_neighbors(True)
	boundary_resis.set_cutoffs(4.9, 2.7)
	boundary_resis.set_dist_exponent(0.7)

	surface_resis = core.select.residue_selector.LayerSelector()
	surface_resis.set_layers(False, False, True) #sets core
	surface_resis.use_sc_neighbors()
	surface_resis.set_use_sc_neighbors(True)
	surface_resis.set_cutoffs(4.9, 2.7)
	surface_resis.set_dist_exponent(0.7)

	core_rlt = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
	core_rlt.aas_to_keep('VILAFM')
	boundary_rlt = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
	boundary_rlt.aas_to_keep('VILAFMWYDERKNQST')
	surface_rlt = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
	surface_rlt.aas_to_keep('DERKNQST')

	layer_design = core.pack.task.operation.DesignRestrictions()
	layer_design.add_selector_rlto_pair(core_resis, core_rlt)
	layer_design.add_selector_rlto_pair(boundary_resis, boundary_rlt)
	layer_design.add_selector_rlto_pair(surface_resis, surface_rlt)

	task_f.push_back(layer_design)
	return(task_f)

#evaluate layer design argparse
false_statements = ['False', 'false', 'FALSE', 'F', 'f', '0', 'Off', 'off']
if args.ld not in false_statements:
	design_layers = True
else:
	design_layers = False

sym = protocols.symmetry.SetupForSymmetryMover()
sym.process_symmdef_file(args.symdef)
sym.apply(pose)

pack = protocols.minimization_packing.symmetry.SymPackRotamersMover(sym_scorefxn)
if design_layers == True:
	pack.task_factory(add_layer_design_to_tf())
else:
	pack.task_factory(set_up_tf())
pack.apply(pose)
pose.dump_pdb(output_name)

