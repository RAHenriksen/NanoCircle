import argparse
import File_parser
import sys
import time

#should you add another command, think Circle-map readextractor, such that its possible to run both Simple and Chimeric
# after each other as to not calculate all reads again. Read which output IP's readextractor returns.

class NanoCircle:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
                usage='''NanoCircle <subprogram> [options]
                Rasmus Amund Henriksen, wql443@alumni.ku.dk, 2020
                Version 1.0.0
                
                The NanoCircle suite
                
                Commands:
                
                Classify     Characterize the different reads into simple or chimeric circle origin 
                Simple       Identifies simple circular DNA, comprised of a single chromosomal fragment
                Chimeric     Identifies chimeric circular DNA, comprised of multiple chromosomal fragments''')

        # create subcommands each with their own arguments
        subparsers = self.parser.add_subparsers()
        self.Classify = subparsers.add_parser(
            name="Classify",
            description='Characterize the different reads into simple or chimeric circle origin',
            prog="NanoCircle Classify",
            usage='''NanoCircle Classify [options]''')

        self.Simple = subparsers.add_parser(
            name="Simple",
            description='Identifies simple circular DNA (single chromosomal fragment)',
            prog="NanoCircle Simple",
            usage='''NanoCircle Simple [options]''')

        self.Chimeric = subparsers.add_parser(
            name="Chimeric",
            description='Identifies chimeric circular DNA (multiple chromosomal fragments)',
            prog="NanoCircle Chimeric",
            usage='''NanoCircle Chimeric [options]''')

        # no commands
        if len(sys.argv) <= 1:
            self.parser.print_help() # prints the help page immediately, but then when using -h option nothing appears
            # Ideally both would work.
        else:
            # Hard coding it so when givin -h it actually prints the hellp message
            if sys.argv[1] == "-h" or sys.argv[1] == "--help":
                self.parser.print_help()

            # The positional arguments
            elif sys.argv[1] == "Classify":
                print("Classify works")
                self.subprogram = self.args_Classify()
                self.args = self.subprogram.parse_args(sys.argv[2:])

            elif sys.argv[1] == "Simple":
                print("Simple works")

                # Defines the subprogram with the arguments given by args_Simple()
                self.subprogram = self.args_Simple()

                # passes all the arguments needed for the simple commands
                self.args = self.subprogram.parse_args(sys.argv[2:])

                #imports the Simple sub-command defined in another script
                import Simple_cmd as Sim
                print("mapq",self.args.mapq)
                # parse the created argument object into the File_parser script to load the data
                # with Input_parser() function
                data_class = File_parser.Input_Parser(self.args.input)
                data_in = data_class.parse_data()

                # passes the loaded data into the Simple_cmd script
                Class_object = Sim.Simple_circ(data_in)
                data_out = Class_object.reverse()

                # passes the results to the File_parser script saving the resulting data
                Out_class = File_parser.Output_Parser(data_out, self.args.output)
                Out_class.dataframe()

            elif sys.argv[1] == "Chimeric":
                print("Chimeric works")
                self.subprogram = self.args_Chimeric()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                import Chimeric_eccDNA as Chim

                Class_object = Chim.Chimeric_circ(self.args.i)
                Class_object.multiply()

            else:
                self.parser.print_help()
                sys.stderr.write(
                    "NanoCircle error: the positional arguments are required\n")
                sys.exit(0)

    def args_Classify(self):
        """
        :return: argument parser for the Simple commands
        """
        parser = self.Classify

        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        # required arguments
        required.add_argument("-o", "--output", required=True, metavar="", help='Tab seperated identified circles')

        # optional arguments
        optional.add_argument("-q", "--mapq", metavar="", default=60, type=int, help='Mapping Quality, default 60')

        # if no arguments are parsed
        if len(sys.argv[2:]) == 0:
            parser.print_help()

        return parser

    def args_Simple(self):
        """
        :return: argument parser for the Simple commands
        """
        parser = self.Simple # The parser defined in __init__

        # creates the argument groups
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        # required arguments
        required.add_argument("-i", "--input",required=True, metavar="", help='Tab seperated potential regions') #hvad med
        required.add_argument("-o", "--output",required=True, metavar="", help='Tab seperated identified circles')

        # optional arguments
        optional.add_argument("-q", "--mapq", metavar="", default=60, type=int, help='Mapping Quality, default 60')

        # if no arguments are parsed
        if len(sys.argv[2:]) == 0:
            parser.print_help()

        return parser

    def args_Chimeric(self):
        """
        :return: argument parser for the Simple commands
        """
        parser = self.Chimeric

        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        # required arguments
        required.add_argument("-t", "--test", required=True, metavar="", help='Tab seperated identified circles')

        # optional arguments
        optional.add_argument("-O", "--OLL", metavar="", default=60, type=int, help='Mapping Quality, default 60')

        # if no arguments are parsed
        if len(sys.argv[2:]) == 0:
            parser.print_help()

        return parser

if __name__ == '__main__':
    NanoCircle()
