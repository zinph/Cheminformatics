#-------------------------------------------------------------------------------
# Name:        MergeFiles
# Author:      zinph
# Created:     29/09/2016
# Copyright:   (c) zinph 2016
#-------------------------------------------------------------------------------
#Objective: To combine multiple files of the same format into one file of either the
#           same format or a different format.
#           Targettig mainly .mol, .sdf, and .txt files.
#How to use it:
#     the program will prompt user for three inputs:
#         (A) file location - you just need to choose the directory in the prompted dialogue       
#         (B) input extension of the files you want to put together (without dot); for example:  mol
#         (C) desired output of the file you want to merge (with extension); for example: newfile.sdf
#-------------------------------------------------------------------------------

from tkinter import filedialog
import os
import os.path

def file_handler():
    '''
    Prompt user to choose a file directory.
    '''
    directory = filedialog.askdirectory()                                       #   ask for directory
    file_list = os.listdir(directory)                                           #   list of all files in the directory assigned to file_list
    i_ext = input("extension of original file (without dot): ")                 #   prompt user for extension of the files s/he wants to compile
    new_file_name = input("name of your merged file (with extension): ")        #   prompt user for the new name file including the desired extension
    os.chdir(directory)                                                         #   change the directory to the selected folder
    if os.path.exists(new_file_name):                                           #   if the file name exists, remove the file.
        os.remove(new_file_name)

    num_files = 0                                                               #   numfiles = variable to keep track of the number of files to be compiled
    for i in file_list:                                                         #   iterate the list of all files
        if i != new_file_name:                                                  #   checkpoint to not encounter file with desired output name
            if i_ext == i[len(i)-len(i_ext):]:                                  #   if extension of file is detected
                num_files += 1                                                  #   increment num_files
                file = open(i,'r')                                              #   open file in read mode
                new_file = open(new_file_name,'a')                              #   open new_file in append mode
                new_file.write(''.join(file.readlines()))                       #   join all the lines from the file and append to new_file
                if 'mol' in new_file_name:
                    new_file.write('$$$$' + '\n' )                                  #   added to comply with .mol format
                new_file.close()
                file.close()                                                    #   close both file and new_file

            elif num_files == 0:
                print ("Files with your requested extension doesn't exist in the folder. ")

    print (str(num_files) + ' of .' + str(i_ext) + ' files are compiled into your desired file with name ' + str(new_file_name))


file_handler()