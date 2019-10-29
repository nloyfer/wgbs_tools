import subprocess
import os
import os.path as op


EXTENSIONS = ('pat.gz', 'pat.gz.csi', 'unq.gz', 'unq.gz.csi', 'beta')


def del_outputs(fname):
    for f in (fname + '.' + ext for ext in EXTENSIONS):
        if op.isfile(f):
            os.remove(f)


def compare_to_ref(fname):
    print('Comparing')
    for f in (fname + '.' + ext for ext in EXTENSIONS):
        cmd = 'diff {f} ref/{f} | head'.format(f=f)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        if p.returncode or error:
            print("diff failed %d %s %s" % (p.returncode, output, error))
            return False
        elif output:
            print('There is a difference in file {}'.format(f))
            return False
    print('all good')
    return True


def compile_ctools():
    r1 = subprocess.call('g++ -std=c++11 ../patter.cpp -o ../patter', shell=True)
    r2 = subprocess.call('g++ -std=c++11 ../match_maker.cpp -o ../match_maker', shell=True)
    return r1 or r2


def run_test(fname, region):
    del_outputs(fname)

    # Compile
    if compile_ctools():
        print('Failed in compilation')
        return

    # Run pipeline
    cmd = 'wgbs_tools bam2pat {}.bam {} --time'.format(fname, region)
    print('Running', cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.returncode != 0:
        print("bam2pat failed %d\n%s\n%s" % (p.returncode, output, error))
        return

    # Compare to reference:
    if not compare_to_ref(fname):
        print('Failed in test for bam2pat\n{}\n{}'.format(output, error))

    # Cleanup
    # del_outputs(fname)


def main():
    print('cd to ', op.dirname(os.path.realpath(__file__)))
    os.chdir(op.dirname(os.path.realpath(__file__)))
    fname = 'test'
    region = ''
    # fname = 'test_chr19'
    # region = '-r chr19'
    run_test(fname, region)
    return


if __name__ == '__main__':
    main()
