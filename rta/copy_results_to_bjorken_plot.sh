

if [ $1 = "generator" ]; then
	echo "Copying results to bjorken_plot/generator"
    cp -r results ../bjorken_plot/generator

elif [ $1 = "ffe_1" ]; then
	echo "Copying results to bjorken_plot/FFE_hydro/xi_0"
    cp -r results ../bjorken_plot/FFE_hydro/xi_0

    elif [ $1 = "ffe_2" ]; then
	echo "Copying results to bjorken_plot/FFE_hydro/xi_9999"
    cp -r results ../bjorken_plot/FFE_hydro/xi_9999

elif [ $1 = "smash" ]; then
	echo "Copying results to bjorken_plot/attractor_new"
    cp -r results ../bjorken_plot/attractor_new
fi


# $1 = generator, ffe_1, ffe_2, attractor