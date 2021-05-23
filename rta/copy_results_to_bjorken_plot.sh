
if [ $1 = "generator" ]; then
	echo "Copying results to bjorken_plot/generator"
    cp -r results ../bjorken_plot/generator

elif [ $1 = "ffe_1" ]; then
	echo "Copying results to bjorken_plot/FFE_hydro/xi_0"
	rm -r ../bjorken_plot/FFE_hydro/xi_0
	mkdir ../bjorken_plot/FFE_hydro/xi_0
    cp -r results ../bjorken_plot/FFE_hydro/xi_0

    elif [ $1 = "ffe_2" ]; then
	echo "Copying results to bjorken_plot/FFE_hydro/xi_9999"
	rm -r ../bjorken_plot/FFE_hydro/xi_9999
	mkdir ../bjorken_plot/FFE_hydro/xi_9999
    cp -r results ../bjorken_plot/FFE_hydro/xi_9999

elif [ $1 = "attractor" ]; then
	echo "Copying results to bjorken_plot/attractor"
    cp -r results ../bjorken_plot/attractor
fi

# $1 = generator, ffe_1, ffe_2, attractor