import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sympy import latex
from sympy import N
from numpy import arange
from essm.variables.units import derive_unit, SI, Quantity
from essm.variables.utils import markdown

def plot_expr2(xvar_min_max, yldata, yllabel=None, yrdata=None,
               yrlabel='', clf=True, npoints=100, ylmin=None, ylmax=None,
               yrmin=None, yrmax=None, xlabel=None,
               colors=None,
               loc_legend_left='best', loc_legend_right='right',
               linestylesl=['-', '--', '-.', ':'], 
               linestylesr=['-', '--', '-.', ':'],
               fontsize=None, fontsize_ticks=None, fontsize_labels=None,
               fontsize_legend=None, xaxispos=None, yaxispos=None,
               fig1=None, legend=True, **args):
    '''
    Plot expressions as function of xvar from xmin to xmax. 
    Expressions to be plotted are given as yldata (for left axis), or yrdata
    (right axis), each of which can be an equation, list of equations or 
    (expr, name) tuples.
    
    
    **Examples:**
    
	from essm.variables import Variable
	from essm.variables.physics.thermodynamics import T_a, nu_a
	from essm.equations.physics.thermodynamics import eq_nua, eq_ka
	vdict = Variable.__defaults__.copy()    
	expr = eq_nua.subs(vdict)
	exprr = eq_ka.subs(vdict)
	xvar = T_a
	yldata = [(expr.rhs, 'full'), (expr.rhs/2, 'half')]
	yrdata = exprr
	# Simplest plot of 1 function
	fig =  plot_expr2((T_a, 273, 373), expr)
	# When plotting 2 functions, need to define label
	fig = plot_expr2((T_a, 273, 373), yldata, yllabel = (nu_a))
	# Plot of two functions on left and one on right axis
	fig = plot_expr2((T_a, 273, 373), yldata, yllabel = (nu_a), yrdata=yrdata)
	# Same as above, but plotting the invers on the right
	fig = plot_expr2((T_a, 273, 373), yldata, yllabel = (nu_a), 
		   yrdata=[Eq(1/exprr.lhs, 1/exprr.rhs)],
		   loc_legend_right='lower right')
	# Modifying fig and ax attributes
	fig.set_dpi(100)
	fig.axes[0].set_ylim(ymin=0)
	fig.axes[0].grid()
    '''
    (xvar, xmin, xmax) = xvar_min_max
    if not colors:
        if yrdata is not None:
            colors = ['black', 'blue', 'red', 'green']
        else:
            colors = ['blue', 'black', 'red', 'green']
    if fontsize:
        fontsize_labels = fontsize
        fontsize_legend = fontsize
        fontsize_ticks = fontsize
    if not fig1:
        plt.close
        plt.clf
        fig = plt.figure(**args)
    else: 
        fig = fig1
    if hasattr(xvar, 'definition'): 
        unit1 = derive_unit(xvar)
        if unit1 != 1:
            strunit = ' (' + markdown(unit1) + ')'
        else: 
            strunit = ''
        if not xlabel:
            xlabel = '$'+latex(xvar)+'$'+ strunit
    else: 
        if not xlabel:
            xlabel = xvar
    if hasattr(yldata, 'lhs'):
        yldata = (yldata.rhs, yldata.lhs)
    if not yllabel:
        if type(yldata) is tuple:
            yllabel = yldata[1]
        else:
            try: 
                yllabel = yldata[0][1]
            except Exception as e1:
                print(e1)
                print('yldata must be equation or list of (expr, name) tuples')
                
    if type(yllabel) is not str: 
        unit1 = derive_unit(yllabel)
        if unit1 != 1:
            strunit = ' (' + markdown(unit1) + ')'
        else: 
            strunit = ''
        
        yllabel = '$'+latex(yllabel)+'$'+ strunit   
    if type (yldata) is not list and type(yldata) is not tuple:
        # If only an expression given
        yldata = [(yldata, '')]
    if type(yldata[0]) is not tuple:
        yldata = [yldata]
    if yrdata is not None:
        if yrlabel == '':
            if hasattr(yrdata, 'lhs'):
                yrlabel = yrdata.lhs 
        if type (yrdata) is not list and type(yrdata) is not tuple:
            # If only an expression given
            yrdata = [yrdata] 
    if type(yrlabel) is not str: 
        yrlabel = '$'+latex(yrlabel)+'$'+ ' (' + markdown(derive_unit(yrlabel)) + ')'            
    
    xstep = (xmax - xmin)/npoints
    xvals = arange(xmin, xmax, xstep)
       
    ax1 =  fig.add_subplot(1, 1, 1)
    if yrdata is not None:
        color = colors[0]
    else:
        color = 'black'
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(yllabel, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    i = 0
    for (expr1, y1var) in yldata:
        linestyle = linestylesl[i]
        if yrdata is None:
            color = colors[i]
        i= i + 1
        try: 
            y1vals = [expr1.subs(xvar, dummy).n() for dummy in xvals]                   
            ax1.plot(xvals, y1vals, color=color, linestyle=linestyle, label=y1var)
        except Exception as e1:
            print([expr1.subs(xvar, dummy) for dummy in xvals])
            print(e1)
    if legend == True and (i > 1 or yrdata is not None):
        plt.legend(loc=loc_legend_left, fontsize=fontsize_legend)
    
    if ylmin is not None:    ax1.set_ylim(bottom=float(ylmin))
    if ylmax is not None:    ax1.set_ylim(top=float(ylmax))
        
    if yrdata is not None:   
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        color = colors[1]
        ax2.set_ylabel(yrlabel, color=color)
        i = 0

        for item in yrdata:
            if type(item) is tuple:   # if item is tuple
                (expr2, y2var) = item
            else:
                try: 
                    (y2var, expr2) = (item.lhs, item.rhs)
                except Exception as e1:
                    print(e1)
                    print('yrdata must be a list of equations or tuples (var, expr)')
                    return
            linestyle = linestylesr[i]
            i = i + 1
            try:
                y2vals = [expr2.subs(xvar, dummy).n() for dummy in xvals]
                ax2.plot(xvals, y2vals, color=color, linestyle=linestyle, label=y2var)
            except Exception as e1:
                print(expr2)
                print([expr2.subs(xvar, dummy).n() for dummy in xvals])
                print(e1)
                
            if not yrlabel:
                if hasattr(yrdata[0], 'lhs'):
                    yrlabel = yrdata[0].lhs

        if type(yrlabel) is not str: 
            yrlabel = '$'+latex(yrlabel)+'$'+ ' (' + markdown(derive_unit(yrlabel)) + ')'       
        ax2.tick_params(axis='y', labelcolor=color)
        if yrmin:    ax2.set_ylim(ymin=float(yrmin))
        if yrmax:    ax2.set_ylim(ymax=float(yrmax))
        leg=ax2.legend(loc=loc_legend_right, fontsize=fontsize_legend)
        ax2.add_artist(leg);
        for item in ([ax2.xaxis.label, ax2.yaxis.label]):
            item.set_fontsize(fontsize_labels)
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

    for item in ([ax1.xaxis.label, ax1.yaxis.label]):
        item.set_fontsize(fontsize_labels)
        
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    if yaxispos is not None:
        print(yaxispos)
        # Make y-axis intersect at xaxispos and remove frame
        ax1.spines['left'].set_position(('data', yaxispos))
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
    if xaxispos is not None:
        print(xaxispos)
        # Make x-axiss intersect at xaxispos and remove frame
        ax1.spines['bottom'].set_position(('data', xaxispos))
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    return fig


