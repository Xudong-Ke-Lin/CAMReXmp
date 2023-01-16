model=$(cd ../../ && echo ${PWD##*/})

l_r=100.0
array=( 8192 4096 2048 1024 512 )
array2=( 0.9 0.9 0.9 0.9 0.9 )

cpu=12

for i in "${!array[@]}"
do
    echo "%s is in %s\n" "${array[i]}" "${array2[i]}"

    mpiexec -n $cpu /data/tycho/xk215/MPHIL_thesis/make_models/$model/./main1d.gnu.MPI.ex inputs l_r=$l_r amr.n_cell=${array[i]}" 10 128" adv.cfl=${array2[i]}  amr.plot_file="/data/tycho/xk215/MPHIL_thesis/"$model"/plt_BrioWu_lr"$l_r"_n_cell"${array[i]}"_cfl"${array2[i]}"_"
done
for f in "/data/tycho/xk215/MPHIL_thesis/"$model"/plt_BrioWu_lr"$l_r"_"*; do ~/AMReX/amrex/Tools/Plotfile/./fextract.gnu.ex $f; done

l_r=1.0
array=( 8192 4096 2048 1024 512 )
array2=( 0.9 0.5 0.5 0.5 0.2 )

for i in "${!array[@]}"
do
    echo "%s is in %s\n" "${array[i]}" "${array2[i]}"

    mpiexec -n $cpu /data/tycho/xk215/MPHIL_thesis/make_models/$model/./main1d.gnu.MPI.ex inputs l_r=$l_r amr.n_cell=${array[i]}" 10 128" adv.cfl=${array2[i]}  amr.plot_file="/data/tycho/xk215/MPHIL_thesis/"$model"/plt_BrioWu_lr"$l_r"_n_cell"${array[i]}"_cfl"${array2[i]}"_"
done
for f in "/data/tycho/xk215/MPHIL_thesis/"$model"/plt_BrioWu_lr"$l_r"_"*; do ~/AMReX/amrex/Tools/Plotfile/./fextract.gnu.ex $f; done

l_r=0.1
array=( 8192 4096 2048 1024 512 )
array2=( 0.01 0.01 0.01 0.005 0.005 )

for i in "${!array[@]}"
do
    echo "%s is in %s\n" "${array[i]}" "${array2[i]}"
    
    mpiexec -n $cpu /data/tycho/xk215/MPHIL_thesis/make_models/$model/./main1d.gnu.MPI.ex inputs l_r=$l_r amr.n_cell=${array[i]}" 10 128" adv.cfl=${array2[i]}  amr.plot_file="/data/tycho/xk215/MPHIL_thesis/"$model"/plt_BrioWu_lr"$l_r"_n_cell"${array[i]}"_cfl"${array2[i]}"_"
done
for f in "/data/tycho/xk215/MPHIL_thesis/"$model"/plt_BrioWu_lr"$l_r"_"*; do ~/AMReX/amrex/Tools/Plotfile/./fextract.gnu.ex $f; done
