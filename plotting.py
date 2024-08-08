import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sympy import latex
from sympy import N
from numpy import arange
from essm.variables.units import derive_unit, SI, Quantity
from essm.variables.utils import markdown

def fun_dict_standardunits(vdict):
    """
    Convert vdict to dictionary consisting of 
    (key, value*standardunit) items, converting
    units, if given, to standardunits, and adding
    standardunits where units were not provided.
    
    vdict is a dictionary with symbolic expressions
    as keys and floats or integers as values, or 
    products of floats or integers with units.
    Returns dictionary with the same keys, but each
    value multiplied by its units.
    """
    for (key, val) in vdict.items():
        if val:
            if hasattr(key, 'definition'):
                targetunit = (key.definition.unit)
            else:
                print("{0} has no defined units.".format(key),
                      "Trying to deduce units from its base SI units",
                      "or assuming it is dimensionless.")
                targetunit = derive_baseunit(key)
            if type(val) == float or type(val) == int:
                val1 = val*targetunit
            else:
                val1 = convert_to(val, targetunit)
                convertedunit = get_unit(val1)
                if convertedunit != targetunit:
                    print("Conversion error for {0}:".format(key),
                          "Target unit was: {0},".format(targetunit),
                          "but final unit was: {0}".format(convertedunit),
                          "Please verify units in original vdict. Aborting.")
                    return vdict
            vdict[key] = val1
    return vdict

def plot_expr2(xvar_min_max, yldata, yllabel=None, yrdata=None,
               yrlabel='', clf=True, npoints=100, ylmin=None, ylmax=None,
               yrmin=None, yrmax=None, xlabel=None, xunit=None, ylunit=None, 
               yrunit=None, colors = ['black', 'red', 'blue', 'green'],
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
    You can provide custom units for axes as xunit and yunit,
    otherwise the default units of xvar and yvar will be assumed. 
    Note that if you provide custom units, you need to substitute values with
    units into expressions before passing them to plot_expr2.
    
    
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
    # Same as above, but plotting the inverse on the right
    fig = plot_expr2((T_a, 273, 373), yldata, yllabel = (nu_a), 
           yrdata=[Eq(1/exprr.lhs, 1/exprr.rhs)],
           loc_legend_right='lower right')
    # Using custom units
    vdict = Variable.__defaults__.copy()
    # Add units to vdict
    vdict = fun_dict_standardunits(vdict)
    vdict[T_a] = 300*kelvin
    expr = eq_Cwa.subs(vdict)
    fig = plot_expr2((P_wa, 0.01, 1), expr)
    fig = plot_expr2((P_wa, 0.01, 1), expr, ylunit= mole/inch**3)
    fig = plot_expr2((P_wa, 0.01, 1), yldata=expr, yrdata=expr, 
                     xunit=bar, ylunit= mole/meter**3, yrunit=mole/inch**3,
                     linestylesr=['--'])
	# Modifying fig and ax attributes
	fig.set_dpi(100)
	fig.axes[0].set_ylim(ymin=0)
	fig.axes[0].grid()
    '''
    def get_unit(expr):
        """Return unit of expression"""
        return expr.subs({x: 1 for x in expr.args if not x.has(Quantity)})

    def get_value(expr):
        """Return value of expression"""
        from sympy import preorder_traversal
        if type(expr) == float or type(expr) == int:
            val = expr
        else:
            units = {arg for arg in preorder_traversal(expr)
                     if isinstance(arg, Quantity)}
            val = expr.subs({x: 1 for x in units})
        return val
    
    (xvar, xmin, xmax) = xvar_min_max

    if fontsize:
        fontsize_labels = fontsize
        fontsize_legend = fontsize
        fontsize_ticks = fontsize

    # Set up figure
    if not fig1:
        plt.close
        plt.clf
        fig = plt.figure(**args)
    else: 
        fig = fig1

    #Set up x-axis label
    try: 
        xunit1 = xunit or derive_unit(xvar)
        if xunit1 != 1:
            strunit = ' (' + markdown(xunit1) + ')'
        else: 
            strunit = ''
        if not xlabel:
            xlabel = '$'+latex(xvar)+'$'+ strunit
    except: 
        if not xlabel:
            xlabel = xvar
            
    # Set up yldata and left label
    try:
        # If Equation provided
        if hasattr(yldata, 'lhs'):
            yldata = [(yldata.rhs, yldata.lhs)]
        # If tuple provided
        if type(yldata) is tuple:
            yldata = [yldata]
        # If only an expression given
        if type (yldata) is not list and type(yldata) is not tuple:             
            yldata = [(yldata, '')]
        for i in range(len(yldata)):
            item = yldata[i]
            # If list of equations given
            if type(item) is not tuple:
                yldata[i] = (item.rhs, item.lhs)        
        if not yllabel:
                yllabel = yldata[0][1]
    except Exception as e1:
        print(e1)
        print('yldata must be equation or list of',
              'equations or (expr, name) tuples')
                
    if type(yllabel) is not str: 
        ylunit1 = ylunit or derive_unit(yllabel)
        if ylunit1 != 1:
            strunit = ' (' + markdown(ylunit1) + ')'
        else: 
            strunit = ''
        
        yllabel = '$'+latex(yllabel)+'$'+ strunit   

    # Set up yrdata and right label
    if yrdata is not None:

        try:
            # If Equation provided
            if hasattr(yrdata, 'lhs'):
                yrdata = [(yrdata.rhs, yrdata.lhs)]
            # If tuple provided
            if type(yrdata) is tuple:
                yrdata = [yrdata]
            # If only an expression given
            if type (yrdata) is not list and type(yrdata) is not tuple:             
                yrdata = [(yrdata, '')]
            for i in range(len(yrdata)):
                item = yrdata[i]
                # If list of equations given
                if type(item) is not tuple:
                    yrdata[i] = (item.rhs, item.lhs)  
            if not yrlabel:
                    yrlabel = yrdata[0][1]
        except Exception as e1:
            print(e1)
            print('yrdata must be equation or list of',
                  'equations or (expr, name) tuples')
                
    if type(yrlabel) is not str: 
        yrunit1 = yrunit or derive_unit(yrlabel)
        if yrunit1 != 1:
            strunit = ' (' + markdown(yrunit1) + ')'
        else: 
            strunit = ''
        
        yrlabel = '$'+latex(yrlabel)+'$'+ strunit       
         
    
   # xdata
    xstep = abs(xmin - xmax)/float(npoints)
    xvals = arange(xmin, xmax, xstep)
    xvals_plot = xvals.copy()  # list of values used for plotting
    if xunit:
        # Convert to standard units for calculation
        baseunit = derive_unit(xvar)
        customunit = xunit
        for i in range(len(xvals)):
            val_orig = xvals[i] * customunit
            val_baseunit = convert_to(val_orig, baseunit)
            xvals[i] = get_value(val_baseunit) # strip units from values

    # Create ax1 with labels
    ax1 =  fig.add_subplot(1, 1, 1)
    if yrdata is not None:
        color = colors[0]
    else:
        color = 'black'
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(yllabel, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    i = 0
    
    # yldata
    for (expr1, y1var) in yldata:
        linestyle = linestylesl[i]
        if yrdata is None:
            color = colors[i]
        i= i + 1
        try:                  
            if ylunit:
                # convert to custom unit
                baseunit = derive_unit(expr1)
                customunit = ylunit
                baseunit_xvar = derive_unit(xvar)
                # Test if units are consistent
                unit1 = get_unit(convert_to(expr1.subs(
                    xvar, xvals[0]*baseunit_xvar),
                    customunit)) 
                if unit1 != customunit:
                    print('WARNING: Requested',
                        'ylunit={0}, but unit of {1} is {2}'.format(
                        ylunit, expr1, unit1))
                
                y1vals = [get_value(convert_to(expr1.subs(
                    xvar, dummy*baseunit_xvar), 
                    customunit)).n() for dummy in xvals]
            else:
                y1vals = [get_value(expr1.subs(xvar, dummy)).n()
                    for dummy in xvals]
            ax1.plot(xvals_plot, y1vals, color=color,
                linestyle=linestyle, label=y1var)
        except Exception as e1:
            print([expr1.subs(xvar, dummy) for dummy in xvals])
            print(e1)
            
            
                
    if legend == True and (i > 1 or yrdata is not None):
        plt.legend(loc=loc_legend_left, fontsize=fontsize_legend)
    
    if ylmin is not None:    ax1.set_ylim(bottom=float(ylmin))
    if ylmax is not None:    ax1.set_ylim(top=float(ylmax))

    # yrdata and ax2
    if yrdata is not None:   
        ax2 = ax1.twinx()  # instantiate a second axes sharing the same x-axis
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
                    print('yrdata must be a list of',
                        'equations or tuples (var, expr)')
                    return
            linestyle = linestylesr[i]
            i = i + 1
            try:
                if yrunit:
                    # convert to custom unit
                    baseunit = derive_unit(expr1)
                    customunit = yrunit
                    baseunit_xvar = derive_unit(xvar)
                    # Test if units are consistent
                    unit1 = get_unit(convert_to(expr2.subs(
                        xvar, xvals[0]*baseunit_xvar), 
                        customunit)) 
                    if unit1 != customunit:
                        print('WARNING: Requested',
                            'yrunit={0}, but unit of {1} is {2}'.format(
                            yrunit, expr1, unit1))
                    
                    y2vals = [get_value(convert_to(
                        expr2.subs(xvar, dummy*baseunit_xvar), 
                        customunit)).n() for dummy in xvals]
                else:
                    y2vals = [get_value(expr2.subs(xvar, dummy)).n()
                        for dummy in xvals]
                ax2.plot(xvals_plot, y2vals, color=color,
                    linestyle=linestyle, label=y2var)
            except Exception as e1:
                print(expr2)
                print([expr2.subs(xvar, dummy).n() for dummy in xvals])
                print(e1)
                
            if not yrlabel:
                if hasattr(yrdata[0], 'lhs'):
                    yrlabel = yrdata[0].lhs

        if type(yrlabel) is not str:
            yrlabel = '${0}$ ({1})'.format(latex(yrlabel),
                markdown(derive_unit(yrlabel)))  
        ax2.tick_params(axis='y', labelcolor=color)
        if yrmin:    ax2.set_ylim(ymin=float(yrmin))
        if yrmax:    ax2.set_ylim(ymax=float(yrmax))
        leg=ax2.legend(loc=loc_legend_right, fontsize=fontsize_legend)
        ax2.add_artist(leg);
        for item in ([ax2.xaxis.label, ax2.yaxis.label]):
            item.set_fontsize(fontsize_labels)
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

    # Tweak font sizes and other details
    for item in ([ax1.xaxis.label, ax1.yaxis.label]):
        item.set_fontsize(fontsize_labels)
        
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    if yaxispos is not None:
        # Make y-axis intersect at xaxispos and remove frame
        ax1.spines['left'].set_position(('data', yaxispos))
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
    if xaxispos is not None:
        # Make x-axiss intersect at xaxispos and remove frame
        ax1.spines['bottom'].set_position(('data', xaxispos))
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    return fig


