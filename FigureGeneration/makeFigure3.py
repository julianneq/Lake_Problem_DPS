from matplotlib import pyplot as plt

def makeFigure3():
    Ax = [1, 2, 3, 5]
    Ay = [7, 6, 5, 4]
    
    Bx = [3, 4, 6]
    By = [7, 2, 1]
    
    light_l_x = [1, 2, 2, 3, 3, 4]
    light_l_y = [7, 7, 6, 6, 5, 5]
    light_u_y = [8, 8, 8, 8, 7, 7]
    
    dark_l_x = [4, 5, 5, 6, 6, 7]
    dark_l_y = [2, 2, 2, 2, 1, 1]
    dark_u_y = [5, 5, 4, 4, 4, 4]
    
    mixed_l_x = [3, 4, 4, 5, 5, 6, 7]
    mixed_l_y = [7, 7, 5, 5, 4, 4, 4]
    mixed_u_y = [8, 8, 8, 8, 8, 8, 8]
    
    ax = plt.subplot(111)
    ax.fill_between(light_l_x, light_l_y, light_u_y, facecolor = '#9ecae1', edgecolor='none')
    ax.fill_between(dark_l_x, dark_l_y, dark_u_y, facecolor = '#4292c6', edgecolor='none')
    ax.fill_between(mixed_l_x, mixed_l_y, mixed_u_y, facecolor = '#084594', edgecolor='none')
    Apts = ax.scatter(Ax, Ay, c='#9ecae1', edgecolor='w', s=100, label='Reference Set A', linewidth=2)
    Bpts = ax.scatter(Bx, By, c='#4292c6', edgecolor='w', s=100, label='Reference Set B', linewidth=2)
    nadir = ax.scatter(7, 8, c='k', edgecolor='w', s=100, label='Reference Point', linewidth=2)
    
    ax.set_xlim(0,7.5)
    ax.set_ylim(0,8.5)
    ax.set_xlabel('Objective 1', fontsize=14)
    ax.set_ylabel('Objective 2', fontsize=14)
    ax.tick_params(axis = 'x', which='both', bottom='off', top = 'off',labelbottom='off')
    ax.tick_params(axis = 'y', which='both', left='off', right = 'off',labelleft='off')
    
    #label hypervolume
    ax.text(2.1, 6.3, 'HV = 5', color='w', fontsize=14, fontweight='bold')
    ax.text(5.1, 3, 'HV = 8', color='w', fontsize=14, fontweight='bold')
    ax.text(5.1, 6.3, 'HV = 12', color='w', fontsize=14, fontweight='bold')
    
    #label ideal point and draw arrows
    ax.scatter(0.5, 0.5, s=500, marker='*', facecolor='#a50f15', edgecolor='k')
    ax.text(0.7, 1.5, 'Ideal', fontsize=14)
    ax.text(0.7, 1, 'Point', fontsize=14)
    ax.arrow(7, 0.5, -6, 0, head_width=0.3, head_length=0.3, fc='k')
    ax.arrow(0.5, 8, 0, -6.75, head_width=0.2, head_length=0.35, fc='k')
    
    #make legend
    A = plt.Rectangle((0, 0), 1, 1, fc='#9ecae1', edgecolor='none')
    B = plt.Rectangle((0, 0), 1, 1, fc='#4292c6', edgecolor='none')
    AB = plt.Rectangle((0, 0), 1, 1, fc='#084594', edgecolor='none')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height*0.85])
    ax.legend([nadir, Apts, Bpts, A, B, AB], ['Reference Point', 'Reference Set A', 'Reference Set B',\
        'Dominated by A only', 'Dominated by B only', 'Dominated by A & B'], scatterpoints = 1, loc='upper center', \
        bbox_to_anchor=(0.5, 1.3), fancybox=True, ncol=2)
    plt.savefig('Figure3.pdf')
    plt.clf()
    
    return None
    
makeFigure3()