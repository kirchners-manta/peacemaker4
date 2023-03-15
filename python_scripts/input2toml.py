''''
This script converts files of the old Peacemaker input format to files 
of toml format.

 - The old Peacemaker input file is a file with the following format:

        # This is a comment
        [section1]
          key1 arg1 arg2
          key2
          key3 arg1 arg2 arg3
        
        [section2]
          # ...
 
 - The toml file has following format:
 
        # This is a comment
        [section1]
          key1 = ["arg1", arg2]
          key2 = []
          key3 = [arg1, "arg2", "arg3"]
 
        [section2]
          # ...

        String arguments are written in double quotes, integer and 
        float arguments are written without quotes.


The script is called from the command line with the following arguments:

        python input2toml.py [old_input_file] [toml_format_file]
'''

import sys

old_input_file = sys.argv[1]
toml_input_file = sys.argv[2]

# open the new toml-format file
with open(toml_input_file, 'w') as f_new:

    # read the old Peacemaker input file
    with open(old_input_file, 'r') as f_old:
        
        # loop over the lines in the old Peacemaker input file
        for line in f_old:

            # remove the newline character from the line
            line = line.rstrip()

            # if the line is empty, write it to the new toml-format file
            if line == '':
                f_new.write(line + '\n')

            # if the line is a comment, write it to the new toml-format file without changing it
            if line.startswith('#'):
                f_new.write(line + '\n')

            # if there is a comment in the end of a line, add it to the new toml-format file
            # without changing it
            elif '#' in line:
                line_without_comment = line.split('#')[0]
                comment = line.split('#')[1]
                                             
                # if the line is a section, write it to the new toml input file
                # without changing it
                if line_without_comment.startswith('['):
                    f_new.write(line_without_comment + ' #' + comment + '\n')

                # if the line is a key, change it to toml format and  write it to the new 
                # toml-format file
                elif line_without_comment.startswith(' '):

                    # extract the key and the arguments
                    key = line_without_comment.split()[0]
                    args = line_without_comment.split()[1:]

                    # convert the arguments to a list, depending on the type of the argument
                    args_list = []
                    for arg in args:
                        try:
                            args_list.append(int(arg))
                        except ValueError:
                            try: 
                                args_list.append(float(arg))
                            except ValueError:
                                args_list.append(str(arg))

                    # if the args_lsit contains only one element, remove the brackets from the list
                    if len(args_list) == 1: 
                        args_list = str(args_list).replace("'", '"').replace('[', '').replace(']', '')

                    # create the toml line
                    toml_line = '  ' + key + ' = ' + str(args_list).replace("'", '"') + ' #' + comment

                    # write the toml line to the new toml-format input file
                    f_new.write(toml_line + '\n')
            
            # if there is no comment in the line and the line is a section, write it  to the new toml-format 
            # file without changing it
            elif line.startswith('['):
                f_new.write(line + '\n')
            
            # if the line is a key, change it to toml format and write it to the new toml-format file
            elif line.startswith(' '):

                # extract the key and the arguments
                key = line.split()[0]
                args = line.split()[1:]

                # convert the arguments to a list, depending on the type of the argument
                args_list = []
                for arg in args:
                    try:
                        args_list.append(int(arg))
                    except ValueError:
                        try: 
                            args_list.append(float(arg))
                        except ValueError:
                            args_list.append(str(arg))

                # if the args_lsit contains only one element, remove the brackets from the list
                if len(args_list) == 1: 
                    args_list = str(args_list).replace("'", '"').replace('[', '').replace(']', '')

                # create the toml line
                toml_line = '  ' + key + ' = ' + str(args_list).replace("'", '"')
                # write the toml line to the new toml-format input file
                f_new.write(toml_line + '\n')

# close the old PM input file
f_old.close()
# close the new toml-format input file
f_new.close()
