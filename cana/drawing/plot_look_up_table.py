import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from IPython.display import display


def plot_look_up_table(n):
    """
    Plot the Look-Up Table of a BooleanNode

    Parameters
    ----------
    n : BooleanNode 
        The BooleanNode to plot the Look-Up Table

    Returns
    -------
    None
    """
    # Check if n.inputs has any values
    if not n.inputs:
        return print('No inputs to plot')
        
    # Init values from BooleanNode
    k = n.k if n.k>=1 else 1
    inputs = n.inputs if not n.constant else [n.name]
    inputlabels = [n.network.get_node_name(i)[0] if n.network is not None else i for i in inputs]
    LUT = n.look_up_table().sort_index(ascending=False)
    # Count number of F in the LUT
    n_fs = LUT.shape[0]
    # Schemata Cell Width and spacing
    cwidth = 60.
    cxspace = 0
    cyspace = 6
    border = 1
    sepcxspace = 21
    # sepcyspace = 15
    dpi = 150.
    # Margins
    top, right, bottom, left, hs = 120, 25, 25, 60, 25
    # Axes Width & Height
    ax1width = ((k*(cwidth+cxspace))+sepcxspace+(cwidth))
    ax1height = (n_fs*(cwidth+cyspace)-cyspace)
    # Figure Width & Height
    fwidth = (left + ax1width + hs + right)
    fheight = (bottom + ax1height + top)
    # Percentages for Axes location
    _ax1w = ((ax1width*100) / fwidth) / 100
    _ax1h = ((ax1height*100) / fheight) / 100
    _bottom = ((bottom*100) / fheight) / 100
    _left = ((left*100) / fwidth) / 100
    _hs = ((hs*100) / fwidth) / 100
    # Init Figure
    fig = plt.figure(figsize=(fwidth/dpi,fheight/dpi), facecolor='w', dpi=dpi)
    ax1 = fig.add_axes((_left,_bottom,_ax1w,_ax1h), aspect=1, label='LUT')

    ### LUT Plot ###

    yticks = []
    patches = []
    x,y = 0.,0.
    #
    for i,r in LUT.iterrows():
        ins = str(r['In:'])
        out = r['Out:']
        x = 0.
        xticks = []
        for input in ins:
            if input == '0':
                facecolor = 'white'
                textcolor = 'black'
            elif input == '1':
                facecolor = 'black'
                textcolor = 'white'      
            text = '{label:s}'.format(label=input)
            ax1.add_artist(Text(x+cwidth/2,y+cwidth/10*4, text=text, color=textcolor, va='center', ha='center',fontsize=14,family='serif'))
            r = Rectangle((x,y), width=cwidth, height=cwidth, facecolor=facecolor, edgecolor='black')
            patches.append(r)
            xticks.append(x+cwidth/2)
            x += cwidth + cxspace

        x += sepcxspace
        is_one = (out == 1) or (out == '1')
        r = Rectangle((x,y), width=cwidth, height=cwidth, facecolor='black' if is_one else 'white', edgecolor='black')
        ax1.add_artist(Text(x-(sepcxspace/2)-(cxspace/2),y+cwidth/10*4, text=':', color='black', va='center', ha='center',fontsize=14,weight='bold',family='serif'))
        ax1.add_artist(Text(x+(cwidth/2),y+cwidth/10*4, text=str(out), color='white' if is_one else 'black', va='center', ha='center',fontsize=14,family='serif'))
        patches.append(r)
        xticks.append(x+cwidth/2)
        yticks.append(y+cwidth/2)
        y += cwidth + cyspace

        #y += sepcyspace

    ax1.add_collection(PatchCollection(patches, match_original=True))
    #
    ax1.set_yticks(yticks)
    ax1.set_yticklabels([r"$f_{%d}$"%(i+1) for i in range(n_fs)[::-1]], fontsize=14)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(inputlabels + ['%s'%(n.name)], rotation=90, fontsize=14)
    #
    ax1.xaxis.tick_top()
    # Remove Tick
    ax1.tick_params(which='major',pad=7)
    for tic in ax1.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    for tic in ax1.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    # Remove Border
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    # Limits
    ax1.set_xlim(-border,ax1width+border)
    ax1.set_ylim(-border,ax1height+border)
    #ax1.invert_yaxis() 

    # FileName
    filename = n.name
    filename = filename.replace('/','_')
    filename = filename.replace(',','_')
    
    ## Display
    display(fig)
    
    plt.close()

def plot_annigen_look_up_table(n):
    """
    Plot annihilation and generation Look-Up Tables side by side.

    For a given BooleanNode `n`, derive two nodes:
    - annihilation: RULE & (NOT X_mid)
    - generation: (NOT RULE) & (X_mid)
    where X_mid is the middle input bit (index k//2), and plot their LUTs.

    Parameters
    ----------
    n : BooleanNode
        The BooleanNode used to derive annihilation and generation LUTs.

    Returns
    -------
    None
    """
    # Guard for nodes without inputs
    if not n.inputs:
        return print('No inputs to plot')

    # Init
    k = n.k if n.k >= 1 else 1
    inputs = n.inputs if not n.constant else [n.name]
    inputlabels = [n.network.get_node_name(i)[0] if n.network is not None else i for i in inputs]

    # Build annihilation and generation nodes from the LUT of n
    lut = n.look_up_table()
    mid_idx = k // 2
    annihilation_outputs_lut = (
        ((lut['Out:'] == '0') & (lut['In:'].str[mid_idx] == '1'))
        .apply(lambda x: '1' if x else '0')
        .tolist()
    )
    generation_outputs_lut = (
        ((lut['Out:'] == '1') & (lut['In:'].str[mid_idx] == '0'))
        .apply(lambda x: '1' if x else '0')
        .tolist()
    )

    from cana.boolean_node import BooleanNode  # local import to avoid cycles at module init
    anni = BooleanNode.from_output_list(annihilation_outputs_lut)
    gen = BooleanNode.from_output_list(generation_outputs_lut)

    LUT_A = anni.look_up_table().sort_index(ascending=False)
    LUT_G = gen.look_up_table().sort_index(ascending=False)

    # Drawing constants (match style of plot_look_up_table)
    n_fs = LUT_A.shape[0]
    cwidth = 60.
    cxspace = 0
    cyspace = 6
    border = 1
    sepcxspace = 21
    dpi = 150.
    # Increase horizontal spacing between subplots to avoid overlap
    # Use at least one cell width plus separator as gap
    top, right, bottom, left = 120, 25, 25, 60
    hs = max(25, int(cwidth + sepcxspace))

    axw = (k * (cwidth + cxspace)) + sepcxspace + (cwidth)
    axh = (n_fs * (cwidth + cyspace) - cyspace)
    fwidth = left + axw + hs + axw + right
    fheight = bottom + axh + top

    _axw = ((axw * 100) / fwidth) / 100
    _axh = ((axh * 100) / fheight) / 100
    _bottom = ((bottom * 100) / fheight) / 100
    _left = ((left * 100) / fwidth) / 100
    _hs = ((hs * 100) / fwidth) / 100

    fig = plt.figure(figsize=(fwidth / dpi, fheight / dpi), facecolor='w', dpi=dpi)
    axA = fig.add_axes((_left, _bottom, _axw, _axh), aspect=1, label='LUT-ANNI')
    axG = fig.add_axes((_left + _axw + _hs, _bottom, _axw, _axh), aspect=1, label='LUT-GEN')

    def _draw_lut(ax, LUT, title):
        yticks = []
        patches = []
        x, y = 0., 0.
        for _, r in LUT.iterrows():
            ins = str(r['In:'])
            out = r['Out:']
            x = 0.
            xticks = []
            for ch in ins:
                if ch == '0':
                    facecolor = 'white'
                    textcolor = 'black'
                else:  # '1'
                    facecolor = 'black'
                    textcolor = 'white'
                ax.add_artist(Text(x + cwidth/2, y + cwidth/10*4, text=ch, color=textcolor,
                                   va='center', ha='center', fontsize=14, family='serif'))
                patches.append(Rectangle((x, y), width=cwidth, height=cwidth, facecolor=facecolor, edgecolor='black'))
                xticks.append(x + cwidth/2)
                x += cwidth + cxspace

            x += sepcxspace
            is_one = (out == 1) or (out == '1')
            patches.append(Rectangle((x, y), width=cwidth, height=cwidth, facecolor='black' if is_one else 'white', edgecolor='black'))
            ax.add_artist(Text(x - (sepcxspace/2) - (cxspace/2), y + cwidth/10*4, text=':', color='black',
                               va='center', ha='center', fontsize=14, weight='bold', family='serif'))
            ax.add_artist(Text(x + (cwidth/2), y + cwidth/10*4, text=str(out), color='white' if is_one else 'black',
                               va='center', ha='center', fontsize=14, family='serif'))
            xticks.append(x + cwidth/2)
            yticks.append(y + cwidth/2)
            y += cwidth + cyspace

        ax.add_collection(PatchCollection(patches, match_original=True))
        ax.set_yticks(yticks)
        ax.set_yticklabels([rf"$f_{{{i+1}}}$" for i in range(n_fs)[::-1]], fontsize=14)
        ax.set_xticks(xticks)
        ax.set_xticklabels(inputlabels + [f'{n.name}'], rotation=90, fontsize=14)
        ax.set_title(str(title), fontsize=12, pad=8)
        ax.xaxis.tick_top()
        ax.tick_params(which='major', pad=7)
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xlim(-border, axw + border)
        ax.set_ylim(-border, axh + border)

    _draw_lut(axA, LUT_A, 'Annihilation')
    _draw_lut(axG, LUT_G, 'Generation')

    display(fig)
    plt.close()

    