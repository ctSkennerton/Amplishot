#!/usr/bin/env python
###############################################################################
#
# Extended PyCogent command line app controller for dealing with positionals
#
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Connor Skennerton"
__copyright__ = "Copyright 2013"
__credits__ = ["Connor Skennerton"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Connor Skennerton"
__email__ = "c.skennerton@gmail.com"
__status__ = "Development"

###############################################################################
import subprocess
import os

from cogent.app.util import CommandLineApplication, ApplicationError
from cogent.app.parameters import Parameter, FlagParameter, ValuedParameter,\
    MixedParameter, Parameters, _find_synonym, is_not_None, FilePath,\
    ParameterError


class RepeatedParameter(Parameter):
    """ A parameter with many occurances on the command line
    """
    def __init__(self, Prefix, Name, Value=None, Delimiter=None, Quote=None,\
        IsPath=False):
        if IsPath and Value:
            Value = map(FilePath, Value)
        super(RepeatedParameter, self).__init__(Name=Name, Prefix=Prefix,\
                Value=Value, Delimiter=Delimiter, Quote=Quote, IsPath=IsPath)
        self._default = Value

    def __str__(self):
        ret = str()
        join_str = self.Id + self.Delimiter + self.Quote
        for i in self.Value:
            ret += ' %s %s%s ' % (join_str, i, self.Quote)
        return ret

    def append(self, value):
        self.Value.append(value)

    def extend(self, value):
        self.Value.extend(value)

    def isOn(self):
        if self.Value is not None:
            return True
        return False

    def isOff(self):
        return not self.isOn()

    def on(self, value):
        if value is not None:
            if isinstance(value, list) or isinstance(value, tuple):
                if self.IsPath:
                    value = map(FilePath, value)
                self.extend(value)
            else:
                if self.IsPath:
                    value = RFilePath(value)
                self.append(value)

        else:
            raise ParameterError,\
            "Turning the ValuedParameter on with value None is the same as "+\
            "turning it off. Use another value."

    def off(self):
        self.Value = None


class ExtentedCommandLneApplication(CommandLineApplication):
    """ Class containing new parameters for dealing with positional arguments
        This class contains a new class variable for positional arguments that
        is a list of Parameters.  There is also a boolean flag which indicates
        whether the positional arguments should be plased at the beginning or
        end of the command line
    """
    
    _positionals = []

    def __init__(self, params=None, positionals=None, InputHandler=None,
            SuppressStderr=None, SuppressStdout=None, WorkingDir=None,
            TmpDir='/tmp', TmpNameLen=20, HALT_EXEC=False,
            PrependPositionals=False):
        """ Set up the ExtendedCommandLineApplication
        """
        super(ExtendedCommandLineApplication, self).__init__(params=params,
                InputHandler=InputHandler, SuppressStderr=SuppressStderr,
                TmpNameLen=TemNameLen, HALT_EXEC=HALT_EXEC)
        self.prependPositionals = PrependPositionals
        if positionals:
            self._postionals = positionals

    def __call__(self, data=None, remove_tmp=True):
        """Run the application with the specified kwargs on data
        
            data: anything that can be cast into a string or written out to
                a file. Usually either a list of things or a single string or 
                number. input_handler will be called on this data before it 
                is passed as part of the command-line argument, so by creating
                your own input handlers you can customize what kind of data
                you want your application to accept

            remove_tmp: if True, removes tmp files
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        out = os.devnull
        err = os.devnull
        if not suppress_stdout:
            out = str(FilePath(self.getTmpFilename(self.TmpDir)))
        if not suppress_stderr:
            err = FilePath(self.getTmpFilename(self.TmpDir))

        outfile = open(out, 'w+b')
        errfile = open(err, 'w+b')

        if data is None:
            input_arg = ''
        else:
            input_arg = getattr(self,input_handler)(data)

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        command = self._command_delimiter.join(filter(None,\
            [self.BaseCommand, str(input_arg)]))
        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command
        # The return value of system is a 16-bit number containing the signal
        # number that killed the process, and then the exit status
        # We only want to keep the exit status so do a right bitwise shift to
        # get rid of the signal number byte
        exit_status = subprocess.call(command, stdout=outfile, stderr=errfile,
                cwd=self.WorkingDir, shell=True)
        #exit_status = system(command) >> 8
        outfile.seek(0)
        errfile.seek(0)
        # Determine if error should be raised due to exit status of
        # appliciation
        if not self._accept_exit_status(exit_status):
            raise ApplicationError, \
             'Unacceptable application exit status: %s\n' % str(exit_status) +\
             'Command:\n%s\nStdOut:\n%s\nStdErr:\n%s\n' % (command,
                                                           outfile.read(),
                                                           errfile.read())


        try:
            result = CommandLineAppResult(\
             outfile, errfile, exit_status,
             result_paths=self._get_result_paths(data))
        except ApplicationError:
            result = self._handle_app_result_build_failure(\
             outfile, errfile, exit_status, self._get_result_paths(data))

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                os.remove(self._input_filename)
                self._input_filename = None

        return result

    def _get_base_command(self):
        """ Returns the full command string 

            input_arg: the argument to the command which represents the input 
                to the program, this will be a string, either 
                representing input or a filename to get input from
         tI"""
        command_parts = []
        # Append a change directory to the beginning of the command to change
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = self.Parameters
        
        command_parts.append(command)
        if self.prependPositionals:
            command_parts.append(self._command_delimiter.join(filter(\
                None, (map(str, parameters.values())))))
            command_parts.extend(map(str, self._positionals))
        else:
            command_parts.extend(map(str, self._positionals))
            command_parts.append(self._command_delimiter.join(filter(\
                None, (map(str, parameters.values())))))

        return self._command_delimiter.join(command_parts).strip()
    
    BaseCommand = property(_get_base_command)
