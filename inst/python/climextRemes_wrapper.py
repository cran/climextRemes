"""
.. module:: climextRemes
   :platform: Unix, Mac, Windows
   :synopsis: Python wrapper for the climextRemes package.

.. moduleauthor:: Christopher Paciorek <paciorek@berkeley.edu>

"""

version = ""
Fort = None

# As of rpy2 >= 3.0.0, use of NULLType rather than RNULLType
import rpy2.rinterface

if hasattr(rpy2.rinterface, 'NULLType'):
    NULLType = rpy2.rinterface.NULLType   # rpy2 >= 3.0.0
else:
    NULLType = rpy2.rinterface.RNULLType  # rpy2 < 3.0.0

def compute_input_map(input_arg_map):
    import rpy2
    import numpy
    import rpy2.robjects.numpy2ri

    if isinstance(input_arg_map, NULLType):
        #result = {}
        #result[name] = None
        return None

    # As of rpy2 >= 3.0.0, return elements of list are numpy arrays not R vectors
    if isinstance(input_arg_map, numpy.ndarray):
        return input_arg_map
    
    # don't change RObject 
    if isinstance(input_arg_map, rpy2.robjects.RObject):
        #result = {}
        #result[name] = input_arg_map
        return input_arg_map

    # named array
    if isinstance(input_arg_map, rpy2.robjects.vectors.ListVector):
        output_names = []
        output_values = []
        try:
            for (i,xname) in enumerate(input_arg_map.names):
                try:
                    value = compute_input_map(input_arg_map[i])
                    output_names.append(xname)
                    output_values.append(value)
                    #print xname, type(input_arg_map[i]), value
                except:
                    output_names.append(xname)
                    output_values.append("Compute Error")
        except:
            output_names.append("Error")
            output_values.append("Error")

        result = dict()
        for (i,f) in enumerate(output_names):
            if isinstance(output_values[i], dict):

                if "__vector" in output_values[i]:
                    result[f + "_names"] = output_values[i]["__vector_names"]
                    result[f] = output_values[i]["__vector"]

                    name_list = result[f + "_names"].tolist()
                    if isinstance(name_list, NULLType):
                        #result[f + "_names"] = None
                        del result[f + "_names"]
                    else:
                        res = result[f + "_names"]
                        for j in range(len(res)):
                            if isinstance(res[j], NULLType):
                                res[j] = None

                elif "__dict" in output_values[i]:

                    if "__dict_row_names" in output_values[i]:
                        result[f + "_row_names"] = output_values[i]["__dict_row_names"]
                        name_list = result[f + "_row_names"].tolist()
                        if isinstance(name_list, NULLType):
                            #result[f + "_row_names"] = None
                            del result[f + "_row_names"]
                        else:
                            res = result[f + "_row_names"]
                            for j in range(len(res)):
                                if isinstance(res[j], NULLType):
                                    res[j] = None
                          
                    if "__dict_column_names" in output_values[i]:
                        result[f + "_column_names"] = output_values[i]["__dict_column_names"]
                        name_list = result[f + "_column_names"].tolist()
                        if isinstance(name_list, NULLType):
                            #result[f + "_column_names"] = None
                            del result[f + "_column_names"]
                        else:
                            res = result[f + "_column_names"]
                            for j in range(len(res)):
                                if isinstance(res[j], NULLType):
                                    res[j] = None

                    if "__dict_dimnames" in output_values[i]:
                        for (i2, f2) in enumerate(output_values[i]["__dict_dimnames"]):
                            result[f + "_dim" + str(i2+1) + "_names"] = f2
                            try:
                                name_list = result[f + "_dim" + str(i2+1) + "_names"].tolist()
                            except:
                                name_list = None
                            if name_list is None or isinstance(name_list, NULLType):
                                #result[f + "_column_names"] = None                                                              
                                del result[f + "_dim" + str(i2+1) + "_names"]
                            else:
                                res = result[f + "_dim" + str(i2+1) + "_names"]
                                for j in range(len(res)):
                                    if isinstance(res[j], NULLType):
                                        res[j] = None
                       
                    result[f] = output_values[i]["__dict"]
                else:
                    result[f] = output_values[i]
            else:
                 result[f] = output_values[i]
        return result

    # put Array and Matrix before Vectors as they also are Vectors and
    # want special handling specific to Array and Matrix
    if isinstance(input_arg_map, rpy2.robjects.vectors.Matrix):
        result = {}
        name = "__dict"
        result[name] = numpy.array(input_arg_map)
        try:
            result[name + "_row_names"] = numpy.array(input_arg_map.rownames)
        except:
            pass
        try:
            result[name + "_column_names"] = numpy.array(input_arg_map.colnames)
        except:
            pass
        return result

    if isinstance(input_arg_map, rpy2.robjects.vectors.Array):
        result = {}
        name = "__dict"
        result[name] = numpy.array(input_arg_map)
        result[name + "_dimnames"] = []
        for (i,f) in enumerate(input_arg_map.dimnames):
            result[name + "_dimnames"].append(compute_input_map(f))
        return result

    # other vectors
    if isinstance(input_arg_map, rpy2.robjects.vectors.IntVector):
        result = {}
        name = "__vector"
        input_arg_map_save = input_arg_map
        try:
            result[name+"_names"] = numpy.array(input_arg_map.names)
        except:
            pass
        has_nan = False
        for f in input_arg_map:
            if f is rpy2.rinterface.NA_Integer:
                has_nan = True
                break
        if has_nan: 
            input_arg_map = [numpy.nan if f is rpy2.robjects.NA_Integer else f for f in input_arg_map]
            result[name] = numpy.array(input_arg_map, numpy.float)
        else:
            result[name] = numpy.array(input_arg_map)
        if name+"_names" not in result.keys() or isinstance(input_arg_map_save.names, NULLType):
            result = result[name]
        return result

    if isinstance(input_arg_map, rpy2.robjects.vectors.FloatVector):
        result = {}
        name = "__vector"
        try:
            result[name+"_names"] = numpy.array(input_arg_map.names)
        except:
            pass
        # Float correctly handles numpy.nan conversion
        result[name] = numpy.array(input_arg_map)
        if name+"_names" not in result.keys() or isinstance(input_arg_map.names, NULLType):
            result = result[name]
        return result

    if isinstance(input_arg_map, rpy2.robjects.vectors.StrVector):
        result = {}
        name = "__vector"
        try:
            result[name+"_names"] = numpy.array(input_arg_map.names)
        except:
            pass
        try:
            result[name] = numpy.array(input_arg_map)
        except:
            result[name] = input_arg_map
        if name+"_names" not in result.keys() or isinstance(input_arg_map.names, NULLType):
            result = result[name]
        return result

    if isinstance(input_arg_map, rpy2.robjects.vectors.BoolVector):
        result = {}
        name = "__vector"
        input_arg_map_save = input_arg_map
        try:
            result[name+"_names"] = numpy.array(input_arg_map.names)
        except:
            pass
        has_nan = False
        for f in input_arg_map:
            if f is rpy2.rinterface.NA_Logical:
                has_nan = True
                break
            if has_nan:
                input_arg_map = [numpy.nan if f is rpy2.robjects.NA_Logical else f for f in input_arg_map]
                result[name] = numpy.array(input_arg_map, numpy.float)
            else:
                result[name] = numpy.array(input_arg_map)
            if name+"_names" not in result.keys() or isinstance(input_arg_map_save.names, NULLType):
                result = result[name]
            return result
    

"""
def compute_input_map(input_arg_map):
    import rpy2
    input_map = None

    if isinstance(input_arg_map, rpy2.rinterface.RNULLType):
        input_map = None
    elif isinstance(input_arg_map, rpy2.robjects.vectors.ListVector):
        input_map = []
        try:
            for (i,name) in enumerate(input_arg_map.names):
                try:
                  value = input_arg_map[i]
                  #print i, name, type(value), isinstance(value, rpy2.robjects.vectors.Vector)
                  if isinstance(value, rpy2.robjects.vectors.ListVector):
                    input_map.append( (name, compute_input_map(value)) )
                  elif isinstance(value, rpy2.robjects.vectors.Vector):
                    input_map.append((name,value[0]))
                  elif isinstance(value, rpy2.rinterface.RNULLType):
                    input_map.append((name,None))
                  else:
                    input_map.append((name,None))
                except:
                  input_map.append((name,None))
        except:
            print "unable to parse", str(input_arg_map)
    elif isinstance(input_arg_map, rpy2.robjects.vectors.Vector):
        input_map = []
        for v in input_arg_map:
          result = compute_input_map(v)
          input_map.append(result)
    elif isinstance(input_arg_map, rpy2.robjects.RObject):
        input_map = input_arg_map # TODO handle this case..
    else: # what else is there?
        input_map = input_arg_map
    
    return input_map

"""

def __initialize_wrapper():
    global version, Fort

    import numpy
    import pandas
    import rpy2
    import rpy2.robjects as robjects
    import rpy2.robjects.help as rh
    from rpy2.robjects import numpy2ri 
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    import fnmatch
    import os
    import re


    # rpy2.robjects.activate()
    if int(rpy2.__version__.split('.')[0]) < 3:
        pandas2ri.activate()
        numpy2ri.activate()
        

    root_dir = os.path.dirname(os.path.realpath(__file__))

    __extRemes__ = importr("extRemes")
    __climextRemes__ = importr("climextRemes")
    __climextRemesHelp__ = rh.Package("climextRemes")

    __version__ = robjects.r('''
    packageVersion("climextRemes")
    ''')[0]

    __version__ = ".".join([str(f) for f in __version__])
    version = __version__

    Fort = robjects.r('''
    data(Fort)
    ord <- order(Fort$year, Fort$month, Fort$day) 
    Fort <- Fort[ord, ]
    ''')

    if int(rpy2.__version__.split('.')[0]) >= 3:
        with rpy2.robjects.conversion.localconverter(rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter):
            Fort = rpy2.robjects.conversion.rpy2py(Fort)
    else:
        Fort = pandas.DataFrame(Fort)
    
    examples = {}
    examples_dir = "{0}/../examples".format(root_dir)
    for file in os.listdir(examples_dir):
        if fnmatch.fnmatch(file, '*_examples.py'):
            filename = examples_dir + "/" + file
            with open(filename, "r") as infile:
                key = file.replace("_examples.py", "")
                value = infile.read()
                examples[key] = value

    arguments_dir = "{0}/../python_help".format(root_dir)
    for file in os.listdir(arguments_dir):
        if fnmatch.fnmatch(file, '*_args.txt'):
            filename = arguments_dir + "/" + file
            with open(filename, "r") as infile:
                key = file.replace("_args.txt", "")
                value = infile.read()
                examples[key + "_arguments"] = value

    #print("examples", examples.keys())

    class climextRemesReturnObject(object):
        def __init__(self):
            pass
        def __repr__(self):
            return str(self.__dict__.keys())

    def __parseResult(obj):
        result = None

        if hasattr(obj, 'names') == False:
            return obj
        
        #print obj.names, rpy2.rinterface.NULL, type(obj), rpy2.rinterface.RNULLType

        if obj.names is rpy2.rinterface.NULL or \
            type(obj.names) is NULLType or \
            (len(obj.names) == 1 and obj.names[0].replace(".","").isdigit()):
            if isinstance(obj, robjects.Vector) and len(obj) == 1:
                result = obj[0]
            else:
                result = robjects.conversion.ri2py(obj)
        elif isinstance(obj, robjects.Vector):
            result = climextRemesReturnObject()
            result.names = numpy.array(obj.names)
            for i in range(len(obj.names)):
                #print obj.names[i]
                result.__dict__[obj.names[i] + "_r"] = numpy.array(obj[i])
                result.__dict__[obj.names[i]] = __parseResult(obj[i])
    
        return result
 
    def __climextRemesWrapFunction(method_name, override_examples):

        names = []
        defaults = []
        signature = robjects.r('''
        library(climextRemes)
        as.list(args({0}))
        '''.format(method_name))

        def quote_args(text, arg_names):
            text = re.sub("( " + " | ".join(arg_names) + " )", "'\\1'", text)
            # As of 0.2.1 changed from curly to straight quote to avoid this error in Py 2.7:
            # SyntaxError: Non-ASCII character '\xe2' in file /accounts/gen/vis/paciorek/R/x86_64-pc-linux-gnu-library/3.5/climextRemes/python/climextRemes_wrapper.py on line 317, but no encoding declared; see http://python.org/dev/peps/pep-0263/ for details
            # This is causing extra spaces because can't tell left vs right single quote, so remove until we go back to curly quotes.
            # text = text.replace("' ", "'")
            # text = text.replace(" '", "'")
            return(text)
        
        for (i,name) in enumerate(signature.names):
            if len(name) > 0:
                name = name.replace(".","_")
                names.append(name)
                try:
                    d = signature[i][0]
                    if isinstance(d, rpy2.robjects.robject.RObject) or isinstance(d, rpy2.rinterface.SexpSymbol):
                        defaults.append(None)
                    else:
                        defaults.append(d)
                except:
                    defaults.append(None)

        gdoc = __climextRemesHelp__.fetch(method_name)
        ## As of rpy2 >= 3.4.0 need \\ in key names
        use_backslash = True if '\\title' in gdoc.sections.keys() else False
        
        key_names = ['title', 'description']
        if use_backslash:
                key_names = ['\\' + nm for nm in key_names]
        gdoc_help = gdoc.to_docstring(key_names)

        if method_name + "_arguments" in override_examples:
            new_arguments = override_examples[method_name + "_arguments"]
            arguments = new_arguments.split("@param ")
            new_arguments = "arguments\n"
            new_arguments += "---------\n"
            argument_names = []
            for arg in arguments[1:]:
                n_arg = arg.split(" ")
                argument_names.append(n_arg[0])
                n_arg[0] = n_arg[0] + ":"
                arg = " ".join(n_arg)
                new_arguments += arg + "\n"

            #new_arguments = new_arguments.replace("@param ", "")
            #gdoc_help += "arguments\n"
            #gdoc_help += "----------\n"
            #gdoc_help += new_arguments #.replace("\n", "\n\n")

            #gdoc_help.replace("\n\ndetails\n", "\n" + new_arguments + "\n\ndetails\n")
            gdoc_help += new_arguments
             
        else:
            key_name = '\\arguments' if use_backslash else 'arguments'
            gdoc_arguments = gdoc.to_docstring([key_name])
            res = gdoc_arguments.split("\n")
            output = ""
            argument_names = []
            for r in res:
                r2 = r.split()
                if len(r2) >= 2:
                    argument_names.append(r2[0])
                    r2[0] = r2[0] + ":"
                    r2 = " ".join(r2)
                    output += " " + r2 + "\n"
                else:
                    output += r + "\n"
            gdoc_help += output + "\n\n"

        #method = getattr(__climextRemes__, method_name)
        #arglist = method.formals().names

        help_keys = list(gdoc.sections.keys())

        key_name = '\\details' if use_backslash else 'details'
        if key_name in help_keys:
            details = gdoc.to_docstring([key_name])
            details = quote_args(details, argument_names)
            gdoc_help += details
            help_keys.remove(key_name)
                
        key_name = '\\value' if use_backslash else 'value'
        if key_name in help_keys:
            value = gdoc.value()
            if use_backslash:    # value needs to be handled differently as of rpy2 >= 3.4.0
                value = ' '.join(["\n\n" if v == "\n" else v.strip() for v in value[1:]])
                value = "\n" + value.strip() + "\n"
            else:
                value = value.replace("\code{", "'")
                value = value.replace("}", "'")        

            gdoc_help += "value\n-----\n\n" + value + "\n"    
            help_keys.remove(key_name)

        if "usage" in help_keys: help_keys.remove('usage')    # this is the R usage so exclude
        if "alias" in help_keys: help_keys.remove('alias')    # this is the R alias so exclude
        if "arguments" in help_keys: help_keys.remove('arguments')    # handled separately
        if "examples" in help_keys: help_keys.remove('examples')    # handled separately
        if "name" in help_keys: help_keys.remove('name')    # this is redundant so exclude

        if "title" in help_keys: help_keys.remove('title')    # handled separately
        if "description" in help_keys: help_keys.remove('description')    # handled separately

        if "\\usage" in help_keys: help_keys.remove('\\usage')    # this is the R usage so exclude
        if "\\alias" in help_keys: help_keys.remove('\\alias')    # this is the R alias so exclude
        if "\\arguments" in help_keys: help_keys.remove('\\arguments')    # handled separately
        if "\\examples" in help_keys: help_keys.remove('\\examples')    # handled separately
        if "\\name" in help_keys: help_keys.remove('\\name')    # this is redundant so exclude

        if "\\title" in help_keys: help_keys.remove('\\title')    # handled separately
        if "\\description" in help_keys: help_keys.remove('\\description')    # handled separately

        others = gdoc.to_docstring(help_keys)
        others = quote_args(others, argument_names)
        gdoc_help += others
        
        if method_name in override_examples:
            gdoc_help += "examples\n"
            gdoc_help += "--------\n"
            gdoc_help += ">>> " + override_examples[method_name].replace("\n", "\n... ")

            #index = gdoc_help.find("examples")
            #if index >= 0:
            #        gdoc_help = gdoc_help[0:index]
            #        if method_name in override_examples:
            #                gdoc_help += "examples\n"
            #                gdoc_help += "---------\n"
            #                gdoc_help += override_examples[method_name].replace("\n", "\n\n")
                
        # some cleanup
        gdoc_help = gdoc_help.replace("title\n-----", "\n", 1)    # 'title' is extraneous
        gdoc_help = gdoc_help.replace("\n \n", "\n\n")
        gdoc_help = gdoc_help.replace("---\n\n\n", "---\n")
        gdoc_help = gdoc_help.replace("\n\n\n", "\n\n")

        help_keys = [hk.replace("\\", "") for hk in help_keys]
        help_keys.extend(["arguments", "examples", "description", "details", "value"])
        for help in help_keys:    # Sphinx thinks repeated characters are section titles
            gdoc_help = gdoc_help.replace(help + "\n" + "-"*len(help), "**" + help    + "**\n")

        gdoc_help = gdoc_help.replace('TRUE', 'True')
        gdoc_help = gdoc_help.replace('FALSE', 'False')
        gdoc_help = gdoc_help.replace('NULL', 'None')

        
        def decorator(*args, **kwargs):
            #print kwargs, args
            new_args = []
            new_kwargs = {}
            for a in args:
                if isinstance(a, dict):
                    new_args.append(rpy2.robjects.vectors.ListVector(a))
                else:
                    if int(rpy2.__version__.split('.')[0]) >= 3 and isinstance(a, numpy.ndarray):
                        # new_args.append(rpy2.robjects.vectors.FloatVector(a)
                        with rpy2.robjects.conversion.localconverter(rpy2.robjects.default_converter + rpy2.robjects.numpy2ri.converter):
                            new_args.append(rpy2.robjects.conversion.py2rpy(a))
                    elif int(rpy2.__version__.split('.')[0]) >= 3 and isinstance(a, pandas.core.frame.DataFrame):
                        with rpy2.robjects.conversion.localconverter(rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter):
                            new_args.append(rpy2.robjects.conversion.py2rpy(a))
                    else:
                        if a is None:
                            a = rpy2.rinterface.NULL
                        new_args.append(a)
            
            for b in kwargs.keys():
                if isinstance(kwargs[b], dict):
                    new_kwargs[b] = rpy2.robjects.vectors.ListVector(kwargs[b])
                else:
                    if int(rpy2.__version__.split('.')[0]) >= 3 and isinstance(kwargs[b], numpy.ndarray):
                        with rpy2.robjects.conversion.localconverter(rpy2.robjects.default_converter + rpy2.robjects.numpy2ri.converter):
                            new_kwargs[b] = rpy2.robjects.conversion.py2rpy(kwargs[b])
                    elif int(rpy2.__version__.split('.')[0]) >= 3 and isinstance(kwargs[b], pandas.core.frame.DataFrame):
                        with rpy2.robjects.conversion.localconverter(rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter):
                            new_kwargs[b] = rpy2.robjects.conversion.py2rpy(kwargs[b])
                    else:
                        if b is None:
                            b = rpy2.rinterface.NULL
                        new_kwargs[b] = kwargs[b]
                
            method = getattr(__climextRemes__, method_name)

            retobj = method(*new_args, **new_kwargs)
        
            #return __parseResult(retobj)
            return compute_input_map(retobj)

        decorator.__doc__ = gdoc_help
        decorator.__name__ = method_name + "_"
        decorator.__names__ = names
        decorator.__defaults__ = tuple(defaults)

        argstr = ", ".join(names)
        fakefunc = "def %s(%s):\n        return %s(%s)\n" % (method_name, argstr, method_name + "_", argstr)
        fakefunc_code = compile(fakefunc, "fakesource", "exec")
        fakeglobals = {}
        eval(fakefunc_code, {method_name + "_": decorator}, fakeglobals)

        f_with_good_sig = fakeglobals[method_name]
        f_with_good_sig.__doc__ = gdoc_help
        f_with_good_sig.__name__ = method_name
        f_with_good_sig.__defaults__ = tuple(defaults)

        return f_with_good_sig

    for __i__ in __climextRemes__._exported_names:
        result = __climextRemesWrapFunction(__i__, examples)
        globals()[__i__] = result

__initialize_wrapper()

