model=$(cd ../../ && echo ${PWD##*/})

cfl=0.2
plot_int=2000
plot_file="/data/melete/xk215/MPHIL_thesis/"$model"/plt_HarrisIM_IM_fully2_cfl"$cfl
mpiexec -n 16 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file"_" num.Maxwell="IM"

ls -1v $plot_file*/Header | tee $plot_file".visit"

