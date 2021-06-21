from pyrosetta import *
from pyrosetta.rosetta import *
import glob

pyrosetta.init('-beta -unmute all -ex1 -ex2aro')
scorefxn = core.scoring.ScoreFunctionFactory.create_score_function('beta')

def set_up_tf(resfile):
	
	limit_aro = protocols.task_operations.LimitAromaChi2Operation()
	limit_aro.chi2max(110)
	limit_aro.chi2min(70)

	tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
	tf.push_back(core.pack.task.operation.InitializeFromCommandline())
	tf.push_back(core.pack.task.operation.IncludeCurrent())
	tf.push_back(limit_aro)
	
	tf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))
	return(tf)

for resfile in glob.glob('*.resfile'):
	design_name = resfile[:4]
	pose = pose_from_file('relax_only_' + design_name + '_0001.pdb')
	#pack = protocols.minimization_packing.symmetry.SymPackRotamersMover(sym_scorefxn)
	pack = protocols.minimization_packing.PackRotamersMover(scorefxn)
	print(resfile)
	pack.task_factory(set_up_tf(resfile))
	pack.apply(pose)
	pose.dump_pdb(resfile.replace('.resfile', '') + '_0001.pdb') 

