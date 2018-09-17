#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {

	if(argc == 1) {
		printf("Digite o nome do dataset!\n");
		return 1;
	}

	else if(argc > 1) {

		for(int c = 1; c < argc;c++) {

			int grid_size[12] = {25, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
			float window_size[5] = {0.05, 0.10, 0.15, 0.20, 0.30};
			char quant[3] = {2, 3, 4};

			char *dataset_name = argv[c];
			char *type = ".txt";

			char *file_name = malloc(strlen(dataset_name)+strlen(type)+1);
			strcat(file_name, dataset_name);
			strcat(file_name, type);

			FILE *file;
			file = fopen(("/home/$USER/experiments/inputs/%s", file_name), "w+");
			
			//AREAFS
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 12; j++) {
					for (int k = 0; k < 5; k++) {
							fprintf(file, "./main areafs fix %d %d %d /home/$USER/experiments/datasets/%s/%s.shp %0.2f\n", grid_size[j], grid_size[j], quant[i], dataset_name, dataset_name, window_size[k]);
					}
				}
			}

			//AREAF
			for(int j = 0; j < 12; j++) {
				for (int k = 0; k < 5; k++) {
					fprintf(file, "./main areaf fix %d %d 1 /home/$USER/experiments/datasets/%s/%s.shp %0.2f\n", grid_size[j], grid_size[j], dataset_name, dataset_name, window_size[k]);
				}
			}

			//MBRC
			for(int j = 0; j < 12; j++) {
				for (int k = 0; k < 5; k++) {
					fprintf(file, "./main mbrc fix %d %d 1 /home/$USER/experiments/datasets/%s/%s.shp %0.2f\n", grid_size[j], grid_size[j], dataset_name, dataset_name, window_size[k]);
				}
			}
		}

	}

	return 0;
}