import os
import csv
import matplotlib.pylab as plt
import numpy as np

def plot_result_hist_size(dataset, query_size):
    are_euler = []
    are_grid = []
    y = []
    hist_size = [10, 20, 30, 50, 70, 100] 

    with open('data.csv','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if (row[4] == 'euler'):
                are_euler.append(float(row[1]))
            elif (row[4] == 'grid'):
                are_grid.append(float(row[1]))
            else:
                print("Nenhum histograma reconhecido")
            y.append(float(row[2]))
    g_label = str(row[4]) + ' histogram'
    plt.plot(hist_size, are_euler, label='euler')
    plt.plot(hist_size, are_grid, label='grid')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Histograma de Euler vs Grade\n\
            dataset=' + str(dataset) + ' query size=' + str(query_size))
    plt.legend() 
    plt.savefig('./graphs/g-' + str(dataset) + '-' + str(query_size) + '.png') 
    csvfile.close()
    plt.clf()


if __name__ == "__main__":
    datasets = ['alertas', 'charminar', 'hidrografia', 'municipios',
            'rodovia']
    hist_size = [10, 20, 30, 50, 70, 100] 
    query_size = [.05, .1, .15, .2, .25, .3]
    for ds in datasets:
         for qs in query_size:
            os.system('rm -rf ./data.csv')
            for hs in hist_size:
                os.system('../main euler areaf fix ' 
                    +str(hs) + ' ' + str(hs) 
                    +' ../../geodatasets/' + ds + '.shp ' + str(qs))
                os.system('../main grid areaf fix ' 
                    +str(hs) + ' ' + str(hs) 
                    +' ../../geodatasets/' + ds + '.shp ' + str(qs))
            plot_result_hist_size(ds, qs)
        
