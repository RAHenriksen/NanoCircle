#!/binf-isilon/common/mfj598/anaconda/bin/python
import argparse
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
                
                Simple       Identifies simple circular DNA, comprised of a single chromosomal fragment
                Chimeric     Identifies chimeric circular DNA, comprised of multiple chromosomal fragments
                Merge        Merge all potential configuration of chimeric eccDNA identified with Chimeric command ''')

        # create subcommands each with their own arguments
        subparsers = self.parser.add_subparsers()

        self.Simple = subparsers.add_parser(
            name="Simple",
            description='Identifies simple circular DNA (single chromosomal fragment)',
            prog="NanoCircle Simple",
            usage='''NanoCircle Simple [options]''')

        self.Chimeric = subparsers.add_parser(
            name="Chimeric",
            description='Identifies configurations of chimeric circular DNA (multiple chromosomal fragments)',
            prog="NanoCircle Chimeric",
            usage='''NanoCircle Chimeric [options]''')

        self.Merge = subparsers.add_parser(
            name="Merge",
            description='Merge all potential configuration of chimeric eccDNA identified with Chimeric command',
            prog="NanoCircle Merge",
            usage='''NanoCircle Merge [options]''')

        # no commands
        if len(sys.argv) <= 1:
            self.parser.print_help() # prints the help page immediately, but then when using -h option nothing appears
            # Ideally both would work.
        else:
            # Hard coding it so when giving -h it actually prints the hellp message
            if sys.argv[1] == "-h" or sys.argv[1] == "--help":
                self.parser.print_help()

            # The positional arguments
            elif sys.argv[1] == "Simple":
                print("Simple works")
                # Defines the subprogram with the arguments given by args_Simple()
                self.subprogram = self.args_Simple()
                # passes all the arguments needed for the simple commands
                self.args = self.subprogram.parse_args(sys.argv[2:])

                #imports the Simple sub-command defined in another script
                import Simple_cmd as Sim

                # passes the loaded data into the Simple_cmd script
                Class_object = Sim.Simple_circ(self.args.input,self.args.ibam,self.args.output,self.args.mapq)
                Class_object.Circle_output()

            elif sys.argv[1] == "Chimeric":
                print("Chimeric works")
                self.subprogram = self.args_Chimeric()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                import Chimeric_cmd as Chim

                Class_object = Chim.Chimeric_circ(self.args.input,self.args.ibam,self.args.output,self.args.mapq)
                Class_object.Region()

            elif sys.argv[1] == "Merge":
                print("Merge works")
                self.subprogram = self.args_Merge()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                import Merge_cmd as Merge

                Class_object = Merge.Merge_config(self.args.input,self.args.output,self.args.dist)
                Class_object.circle_df()

            else:
                self.parser.print_help()
                sys.stderr.write(
                    "NanoCircle error: the positional arguments are required\n")
                sys.exit(0)

    def args_Simple(self):
        """
        :return: argument parser for the Simple commands
        """
        parser = self.Simple # The parser defined in __init__

        # creates the argument groups
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        # required arguments
        required.add_argument("-i", "--input",required=True, metavar="", help='Tab seperated potential regions')
        required.add_argument("-b", "--ibam", required=True, metavar="", help='bamfile')
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
        required.add_argument("-i", "--input",required=True, metavar="", help='Tab seperated potential regions')
        required.add_argument("-b", "--ibam", required=True, metavar="", help='bamfile')
        required.add_argument("-o", "--output", required=True, metavar="",
                              help='Tab seperated identified circles, No. columns corresponds to most complex chimeric circle')
        # optional arguments
        optional.add_argument("-q", "--mapq", metavar="", default=60, type=int, help='Mapping Quality, default 60')

        # if no arguments are parsed
        if len(sys.argv[2:]) == 0:
            parser.print_help()

        return parser

    def args_Merge(self):
        """
        :return: argument parser for the Simple commands
        """
        parser = self.Merge

        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        # required arguments
        required.add_argument("-i", "--input", required=True, metavar="",
                              help='Tab seperated identified chimeric eccDNA configurations')
        required.add_argument("-o", "--output", required=True, metavar="",
                              help='Tab seperated merged chimeric eccDNA configurations')

        # optional arguments
        optional.add_argument("-d", "--dist", metavar="", default=1000, type=int, help='Distance between potential overlapping regions')

        # if no arguments are parsed
        if len(sys.argv[2:]) == 0:
            parser.print_help()

        return parser

if __name__ == '__main__':
    NanoCircle()
