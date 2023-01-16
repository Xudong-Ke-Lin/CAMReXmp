model=$(cd ../../ && echo ${PWD##*/})

l_r=1.0
cfl=0.9
n_cell=512
plotfile="/data/melete/xk215/MPHIL_thesis/"$model"/plt_BrioWu2_hyp_im_MUSTA1SLIC_lr"$l_r"_n_cell"$n_cell"_cfl"$cfl"_"
mpiexec -n 12 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp1d.gnu.MPI.ex inputs l_r=$l_r amr.n_cell=$n_cell" 10 128" adv.cfl=$cfl amr.plot_file=$plotfile

for f in $plotfile*/; do ~/AMReX/amrex/Tools/Plotfile/./fextract.gnu.ex $f; done

