model=$(cd ../../ && echo ${PWD##*/})

cfl=0.1
plot_int=10000
#plot_int=0
plot_file="/data/melete/xk215/MPHIL_thesis/"$model"/plt_Harris_WENO_divFree_RK3_Strang2_correctBC_manualBC_veryLowRes_cfl"$cfl
mpiexec -n 16 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file"_" amr.n_cell="128 64 0"  #stop_time=25.6037796
ls -1v $plot_file*/Header | tee $plot_file".visit"
