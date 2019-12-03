import os

target_location = 'autoFiles/'
result_file_name = 'results.txt'
summary_file_name = 'results_summary.txt'
base_python_command = 'python {LOCATION}/baseFiles/{TYPE}_script.py {PW} {KPT}'
base_submission_file = 'baseFiles/script.sh'
submission_file_name = 'script.sh'

max_submissions = 999 # reduce for testing, increase to run as many as possible.
pw_range = range(100, 400, 50)
k_range = range(2, 8)

def generate_filetree():
    os.mkdir(target_location)
    os.chdir(target_location)

    for i in pw_range:
        pw_locaiton = 'pw{0}/'.format(i)
        os.mkdir(pw_locaiton)
        os.chdir(pw_locaiton)

        for j in k_range:
            k_location = 'k{0}/'.format(j)
            os.mkdir(k_location)
            os.chdir(k_location)

            for k in ['co', 'slab', 'comb']:
                os.mkdir(k)
                os.chdir(k)

                f = open('stdout', 'w')
                f.write('job not yet run')
                f.close()

                os.chdir('..')
            os.chdir('..')
        os.chdir('..')
    os.chdir('..')


def generate_scripts():
    input = open(base_submission_file, 'r')
    script_text = input.read()
    input.close()
    base_location = os.getcwd()

    os.chdir(target_location)

    for i in pw_range:
        os.chdir('pw{0}/'.format(i))

        for j in k_range:
            os.chdir('k{0}/'.format(j))
            for k in ['co', 'slab', 'comb']:
                os.chdir('{0}/'.format(k))

                job_name = 'pw{0}k{1}{2}'.format(i, j, k)
                sh_script_location = os.getcwd()
                python_command = base_python_command.format(LOCATION = base_loc$

                output = open(submission_file_name, 'w')
                sh_script_location = os.getcwd()
                python_command = base_python_command.format(LOCATION = base_loc$

                output = open(submission_file_name, 'w')
                output.write(script_text.format(NAME = job_name,
                    SCRIPT_LOCATION = sh_script_location,
                    COMMAND = python_command))
                output.close()

                os.chdir('..')
            os.chdir('..')
        os.chdir('..')
    os.chdir('..')

def run_scripts():
    submissions = 0
    os.chdir(target_location)

    for i in pw_range:
        os.chdir('pw{0}/'.format(i))

        for j in k_range:
            os.chdir('k{0}/'.format(j))

            for k in ['co', 'slab', 'comb']:
                os.chdir('{0}/'.format(k))

                if not os.path.isfile(result_file_name) and os.path.isfile('std$
                    os.remove('stdout')
                    print('Missing: pw{0}k{1}{2}'.format(i, j, k))
                    os.system('qsub ' + submission_file_name)
                    submissions += 1

                os.chdir('..')
            os.chdir('..')
        os.chdir('..')
    os.chdir('..')

def gather_results():
    none_missing = True
    all_results = ''

    os.chdir(target_location)

    for i in pw_range:
        os.chdir('pw{0}/'.format(i))

        for j in k_range:
            os.chdir('k{0}/'.format(j))

            group_results = []
            for k in ['co', 'slab', 'comb']:
                os.chdir('{0}/'.format(k))

                if os.path.isfile(result_file_name):
                    r = open(result_file_name, 'r')
                    group_results.append(float(r.read()))
                    r.close()
                else:
                    none_missing = False

                os.chdir('..')

            if len(group_results) is 3:
                E_ads = group_results[2] - group_results[0] - group_results[1]
                all_results += 'pw{0}k{1}: {2}\n'.format(i, j, E_ads)

            os.chdir('..')
        os.chdir('..')
    os.chdir('..')

    output = open(summary_file_name, 'w')
    output.write(all_results)
    output.close()
    if none_missing:
        print('All results complete.')

def main():
    if not os.path.isdir(target_location):
        generate_filetree()
        print('Directory generated - rerun to submit jobs.')
    else:
    generate_scripts()
        run_scripts()
        gather_results()


if __name__== '__main__':
  main()
