import sys
from enum import Enum

# required for python 3.7 or higher
if sys.version_info < (3, 7):
    print("ERROR: Python 3.7 or higher is required.")
    sys.exit(1)
import time
import argparse
import datetime
import getpass
import json
import os
import pathlib
import re
import subprocess

debug = False

VERSION = "1.0.9"

MSK_ACCESS_ANALYSIS_TYPE_IDS = [611407000]


class Sequencer(Enum):
    OTHER = 6000
    ROCHE_FLX = 106000
    ROCHE_JUNIOR = 206000
    ILLUMINA_HISEQ1000 = 306000
    ILLUMINA_HISEQ = 406000
    ILLUMINA_GAIIX = 506000
    ILLUMINA_MISEQ = 606000
    ILLUMINA_HISCANSQ = 706000
    LIFE_TECH_5500 = 806000
    LIFE_TECH_5500XL = 906000
    LIFE_TECH_IONTORRENT = 1006000
    LIFE_TECH_IONPROTON = 1106000
    ILLUMINA_NEXTSEQ = 1306000
    ILLUMINA_MISEQSINGLEENDED = 1406000
    ILLUMINA_MINISEQ = 1506000
    ILLUMINA_NOVASEQ = 1606000
    LIFE_TECH_IONS5 = 1706000
    MGI_DNBSEQ_G400 = 1806000
    MGI_DNBSEQ_G50 = 1906000
    MICROARRAY = 100106000
    UNKNOWN = 000000000

    @staticmethod
    def from_value(value: int):
        for sequencer in Sequencer:
            if sequencer.value == value:
                return sequencer
        print(f"Error: Unknown sequencer ID {value}")
        return sequencer.UNKNOWN

    # short name for the sequencer
    def short_name(self):
        return self.name.split(".")[-1]


def is_msk_pipeline(_analysis_type_id):
    return _analysis_type_id in MSK_ACCESS_ANALYSIS_TYPE_IDS


class PatientHelper:
    ONE_PATIENT_RE = r"{\"medicalInformationId\":\d+,\"personalInformationId\":\d+,\"userRef\":\".+\"}"
    ALL_PATIENTS_RE = rf"(\[{ONE_PATIENT_RE}(,{ONE_PATIENT_RE})*\])"

    MAX_P_REF_LENGTH = 70

    @staticmethod
    def read(path, recurse, regex_file):
        """
        :param path: Folder of fastq files
        :param recurse:
        :param regex_file:
        :return: dict of {patient_ref_01: {sample_id_01: [tag, file1, file2], sample_id_02: [tag, file1, file2]}, . . .}
        """

        print("Reading FastQ folder")

        if os.path.exists(path):
            patient_data = PatientHelper._sort_patients(PatientHelper._read(path, recurse, regex_file))
        else:
            print(f"Error: Couldn't find {path}")
            patient_data = None

        return patient_data

    @staticmethod
    def _read(path, recurse, regex_file):
        """
        Get the patient info of fastq files in path. If recurse is True, recurse to subfolders
        :param path:
        :param recurse:
        :param regex_file:
        :return:
        """

        # Check sub folders if recursion is set
        if recurse:
            fastq_files = list(pathlib.Path(path).rglob("*"))
        else:
            fastq_files = list(pathlib.Path(path).glob("*"))

        if len(fastq_files) == 0:
            print("Error: no files found")
            patient_data = None
        else:
            print(f"Found {len(fastq_files)} files")
            patient_data = PatientHelper._get_patient_data(fastq_files, regex_file)

        return patient_data

    @staticmethod
    def _get_patient_data(filenames, regex_override=None):
        """
        :param filenames: list of fastq filenames
        :return: dict of {patient_ref: sample_id1: [tag, file1, file2], sample_id2: [tag, file1, file2], . . . }
        """

        patients = dict()
        regex_helper = RegexHelper(regex_override)
        for file in filenames:
            file = str(file)
            updated = regex_helper.update(file, patients)
            if not updated:
                print(f"Skipping {file}")

        return patients

    @staticmethod
    def update_patient(patients, p_ref, s_id, tag, file):
        """
        Add/update patient info in patients dict
        :param patients:
        :param p_ref:
        :param s_id:
        :param file:
        :param tag:
        :return:
        """
        if p_ref in patients:
            if s_id in patients[p_ref]:
                patients[p_ref][s_id].append(file)
            else:
                patients[p_ref][s_id] = [tag, file]
        else:
            patients[p_ref] = {s_id: [tag, file]}

    @staticmethod
    def _sort_patients(patients):
        """
        Return a copy of patients with the samples in sorted order
        :param patients:
        :return:
        """

        if patients is None:
            return None

        ids2p_ref = dict()
        for p_ref, samples in patients.items():
            for _id in samples.keys():
                ids2p_ref[_id] = p_ref
        sorted_ids = list(ids2p_ref.keys())
        sorted_ids.sort(key=lambda x: int(x[1:]))

        sorted_patients = dict()
        for _id in sorted_ids:
            p_ref = ids2p_ref[_id]
            if p_ref in sorted_patients:
                sorted_patients[p_ref][_id] = patients[p_ref][_id]
            else:
                sorted_patients[p_ref] = {_id: patients[p_ref][_id]}

        return sorted_patients

    @staticmethod
    def create_patients(command, patient_data, client_id, forceplatform):
        """
        :param command:
        :param patient_data:
        :param client_id:
        :return: dict of {patient_ref: (personalInformationId, medicalInformationId), . . . }
        """

        patient_refs = "--patient-ref=" + ",".join(patient_data.keys())

        # Build the create command
        print("Creating patients")
        if forceplatform:
            create_command = command[:] + ["patient", "-c", patient_refs, "-fp"]
        else:
            create_command = command[:] + ["patient", "-c", patient_refs]

        if client_id:
            _ = AuthHelper.access_alt_client(create_command, client_id)
        else:
            result_create = run(create_command)
            _ = result_create.stdout
            Logger.log(result_create)

        # Build the list command
        print("Listing patients")
        if forceplatform:
            list_command = command[:] + ["patient", "-l", patient_refs, "-fp"]
        else:
            list_command = command[:] + ["patient", "-l", patient_refs]

        if client_id:
            stdout_list = AuthHelper.access_alt_client(list_command, client_id)
        else:
            result_list = run(list_command)
            stdout_list = result_list.stdout
            Logger.log(result_list)

        return PatientHelper._extract_patients(stdout_list)

    @staticmethod
    def _extract_patients(stdout):
        # dict of IDs to return
        ids = None

        # Find the returned patient list and convert to Json
        patient_list = re.search(PatientHelper.ALL_PATIENTS_RE, stdout)
        if patient_list:
            try:
                print("Parsing IDs")
                result_json = json.loads(patient_list.groups()[0])
                ids = dict()
                for patient in result_json:
                    patient_ref = patient['userRef']
                    per_id = patient['personalInformationId']
                    med_id = patient['medicalInformationId']
                    ids[patient_ref] = (per_id, med_id)
            except json.decoder.JSONDecodeError:
                print("Error: couldn't decode patient IDs")
        else:
            print("Error: couldn't find patient list in stdout")

        return ids

    @staticmethod
    def is_valid(patient_data):
        """
        Checks that the patient data conforms to the expected format
        :param patient_data:
        :return:
        """

        print("Validating patient data ...")

        def _is_mys(_tags):
            return "D" in _tags and "R" in _tags

        def _is_tumornormal(_tags):
            return "N" in _tags and "T" in _tags

        def _is_multilibrary(_tags):
            return "lib1" in _tags and "lib2" in _tags

        def _is_msk(_tags):
            return "CP" in _tags or "CN" in _tags

        _valid = True

        for p_ref, samples in patient_data.items():
            probs = []

            # Check length of p_ref
            if len(p_ref) > PatientHelper.MAX_P_REF_LENGTH:
                probs.append(f"Patient refs can be no longer than {PatientHelper.MAX_P_REF_LENGTH} characters")

            tags = [samples[s_id][0] for s_id in samples.keys()]
            num_samples = len(samples)
            if num_samples == 2:
                # Should be DNA/RNA, Normal/Tumour, or unspecified
                if not _is_mys(tags) and not _is_tumornormal(tags) and not _is_multilibrary(tags) and not _is_msk(
                        tags) and tags != ["", ""]:
                    probs.append("Patients with two samples should have D and R (mys), T and N (tumorNormal), or none")
                    probs[-1] += f" - found <{tags[0]}> and <{tags[1]}>"
            elif num_samples == 1:
                # One sample should have D, R, or no tag
                if tags[0] not in ["", "D", "R", "T", "N", "CN", "CP", "lib1", "lib2"]:
                    probs.append(f"Single sample patients should be D, R, or empty - found {tags[0]}")
            else:
                probs.append(f"Each patient should have one or two samples - found {num_samples}")

            # If errors were found, inform the user
            if len(probs) > 0:
                _valid = False
                print(f"Error{'s' if len(probs) > 1 else ''} for patient ref {p_ref}")
                for p in probs:
                    print(p)
                print(json.dumps(patient_data[p_ref], indent=3))

        return _valid

    @staticmethod
    def confirm_normal_sample_analysis(patient_data):
        """
        Check if there's a normal sample in the patient data and ask for user's confirmation
        :param patient_data: patient data
        :return: True if user confirms or there's no normal sample, False otherwise
        """
        for samples in patient_data.values():
            for data in samples.values():
                if data[0] == 'N':
                    # Display the message and ask for confirmation
                    confirmation = input(
                        "I acknowledge that when a Normal sample is analyzed, incidental findings of germline variants may be reported. Do you wish to continue? (y/n) ").strip()[
                        0].lower()
                    if confirmation != 'y':
                        print("User didn't confirm. Stopping the execution.")
                        return False
                    return True
        return True

    @staticmethod
    def is_valid_msk(patient_data):
        """
        Checks that the patient data conforms to the expected format
        :param patient_data:
        :return:
        """

        print("Validating MSK patient data ...")

        _valid = True
        tumor_sample_exists = False

        for p_ref, samples in patient_data.items():
            probs = []
            tags = [samples[s_id][0] for s_id in samples.keys()]

            # Check all samples must have a -CP, -CN -T or -N naming convention
            if any(tag not in ['T', 'CP', 'CN', 'N'] for tag in tags):
                probs.append(
                    'MSK Products require all samples to contain one of the following suffixes: -T, -CP, -CN or -N.')

            # Check if a -T sample exists
            if 'T' in tags:
                tumor_sample_exists = True

            # Check unpaired samples can be -T, -CP, -CN but not -N
            if 'N' in tags and 'T' not in tags:
                probs.append('MSK Products do not allow uploading unmatched normal samples (-N without respective -T)')

            # If errors were found, inform the user
            if len(probs) > 0:
                _valid = False
                print(f"Error{'s' if len(probs) > 1 else ''} for patient ref {p_ref}")
                for p in probs:
                    print(p)
                print(json.dumps(patient_data[p_ref], indent=3))

        if not tumor_sample_exists:
            _valid = False
            print('MSK Products require at least one tumor sample (-T) across all samples')

        print(f"MSK sample name validation complete. {'Valid' if _valid else 'Invalid'}")
        return _valid


class AuthHelper:

    @staticmethod
    def access_alt_client(command, client_id):
        """
        Authorise secondary client for patient command
        :param command:
        :param client_id:
        :return:
        """
        stdout = ""

        # Ask user for username/password and add to command
        user, password = AuthHelper._get_auth_info()
        command += ["--client-id", client_id, "-u", user, "-p", password]

        # Start uploader
        with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE) as auth_process:
            # Try to read in token prompt
            prompt = AuthHelper.get_prompt(auth_process)
            if AuthHelper.contains_coordinates(prompt):
                token = input(prompt)
                print("Sending token to subprocess")
                stdout, stderr = auth_process.communicate(token.encode())
                stdout = stdout.decode()
                if stderr:
                    print(f"stderr: {stderr}")

        return stdout

    @staticmethod
    def get_prompt(auth_process):
        """
        Read the prompt one byte at a time until it finds the colon character.
        Assumes the terminal colon is the only one. Current prompt is "Please enter token for coordinates [1, A]: "
        :param auth_process:
        :return:
        """

        prompt = ""
        bytes_to_read = 1

        next_char = auth_process.stdout.read(bytes_to_read).decode()
        if len(next_char) == bytes_to_read:
            while next_char != ":":
                prompt += next_char
                next_char = auth_process.stdout.read(bytes_to_read).decode()
        else:
            print("Error when reading prompt from subprocess")

        return f"{prompt}: "

    @staticmethod
    def _get_auth_info():
        """
        Ask user for username and password
        :return:
        """

        user = input("Please enter your username: ")
        password = getpass.getpass("Please enter your password: ")

        return user, password

    @staticmethod
    def contains_coordinates(prompt):
        """
        Return whether the line contains coordinates (e.g. [1, A]) to confirm we have the token prompt
        :param prompt:
        :return:
        """
        pattern = r"\[[1-8], [A-H]\]"
        match = re.search(pattern, prompt)
        return match


class UserHelper:
    @staticmethod
    def get_user_info(command):
        """
        :return: userId and clientID from userInfo command
        """

        global debug
        user_id = -1
        client_id = -1

        # Execute `userInfo`
        command += ["userInfo"]
        result = run(command)
        if debug:
            print("Fetching userInfo")
            Logger.log(result)
        try:
            result_json = json.loads(result.stdout)
            if "userId" in result_json:
                user_id = result_json["userId"]
            if "clientId" in result_json:
                client_id = result_json["clientId"]
        except json.decoder.JSONDecodeError:
            print("Error: couldn't retrieve userInfo")

        return user_id, client_id

    @staticmethod
    def get_pipeline(command, pipeline_id, forceplatform):
        """
        If no pipeline_id provided, list those found and user chooses. If valid choice return associated sequencerId
        :return: pipelineId and sequencerId for provided pipeline_id (default to -1, -1)
        """
        sequencer_id = -1
        analysis_type = ""

        # Execute the `pipeline --list` command and get the output as Json
        print("Fetching available pipelines")
        if forceplatform:
            command += ["pipeline", "--list", "-fp"]
        else:
            command += ["pipeline", "--list"]
        result = run(command)
        result_json = json.loads(result.stdout)

        # A single pipeline found; put it in a list
        if isinstance(result_json, dict):
            result_json = [result_json]

        # Create map of pipeline_id -> (pipeline_name, sequencer_id)
        available_pipelines = dict()
        for pipeline in result_json:
            if debug:
                print(pipeline)
            pid = pipeline.get("pipeline_id", pipeline.get("id"))
            name = pipeline.get("pipeline_name", pipeline.get("dgName"))
            seq = pipeline.get("sequencer_id", pipeline.get("sequencerId"))
            at = pipeline.get("analysis_type_id", pipeline.get("analysisTypeId"))
            available_pipelines[pid] = (name, seq, at)

        # Pipeline not specified on command line; ask user now
        if pipeline_id == -1:
            if len(available_pipelines) == 0:
                print("No pipelines found")
                exit(3)
            print("Multiple pipelines available")
            for pid, (name, seq, _) in available_pipelines.items():
                print(f'{pid:5}: {name} ({Sequencer.from_value(seq).short_name()})')
            pipeline_id = input("Enter pipeline ID: ")

        # Convert the pipeline_id to int
        if pipeline_id.isnumeric():
            pipeline_id = int(pipeline_id)

        # Check available pipelines for specified id
        for pid, (_, seq, at) in available_pipelines.items():
            if pipeline_id == pid:
                sequencer_id = seq
                analysis_type = at  # This is the analysis type of the pipeline
                break

        if sequencer_id == -1:
            print(f"Error: pipeline {pipeline_id} is invalid")

        return pipeline_id, sequencer_id, analysis_type


class JsonBuilder:
    # Map tags to analyses->definition->libraryType
    lib_type = {
        "R": "rna",
        "D": "dna",
        "N": "",
        "T": "",
        "": "dna",
        "CP": "",
        "CN": "",
        "lib1": "",
        "lib2": ""
    }

    # Map tags to topology->references->role
    role = {
        "D": "dna",
        "R": "rna",
        "N": "normal",
        "T": "tumor",
        "lib1": "lib1",
        "lib2": "lib2"

    }

    # Map tags to topology->definition->type
    # (R and T _should be_ redundant as the file names are sorted)
    top_type = {
        "D": "mys",
        "R": "mys",
        "N": "tumorNormal",
        "T": "tumorNormal",
        "lib1": "multiLibrary",
        "lib2": "multiLibrary"
    }

    # Tag to sample type id
    # If the pipeline selected is MSK-ACCESS the sample type is:
    # Tumor (suffix -T) > default to cfDNA 1308000
    # Normal (suffix -N) > default to WBC DNA 1408000
    # Controls (suffix -CP or -CN) > default to N/A 1508000

    tag_to_sample_type_id_map = {
        "R": 108000,  # RNA
        "D": 108000,  # DNA
        "N": 108000,  # Normal
        "T": 108000,  # Tumor
        "lib1": 108000,  # DNA
        "lib2": 108000  # DNA
    }

    tag_to_sample_type_id_map_msk_access = {
        "R": 108000,  # RNA
        "D": 108000,  # DNA
        "N": 1408000,  # WBC DNA
        "T": 1308000,  # cfDNA
        "CP": 1508000,  # N/A
        "CN": 1508000  # N/A
    }

    @staticmethod
    def build_json(user_ref, user_id, client_id, pipeline_id, sequencer_id, patient_data, patient_ids, sample_type_id,
                   _is_msk_pipeline, bds_mapping):
        """
        :param user_ref:
        :param user_id:
        :param client_id:
        :param pipeline_id:
        :param sequencer_id:
        :param patient_data:
        :param patient_ids:
        :param sample_type_id:
        :return: ADE format json object
        """
        ade_dict = {
            "protocolName": "ADE",
            "protocolVersion": "1",
            "client": {
                "id": client_id,
                "userId": user_id
            },
            "request": {
                "definition": {
                    "userRef": user_ref,
                    "sequencerId": sequencer_id,
                    "requestDate": int(time.time()),
                    "isPairedEnd": True,
                    "isPrevent": False
                },
                "state": None,
                "analyses": JsonBuilder._build_analyses(patient_data, patient_ids, pipeline_id, sample_type_id,
                                                        _is_msk_pipeline, bds_mapping),
                "topology": JsonBuilder._build_topologies(patient_data),
                "files": []
            }
        }

        return json.dumps(ade_dict, indent=3)

    @staticmethod
    def _build_analyses(patient_data, patient_ids, pipeline_id, sample_type_id, _is_msk_pipeline, _bds_mapping):
        """
        :param patient_data:
        :param pipeline_id:
        :return: list of analyses for ADE
        """
        analyses = []
        # order analysis by MID numerical value e.g. S1, S2, S3
        # Flatten the patient data, so we have one level of samples per key
        flattened_data = []
        for patient, samples in patient_data.items():
            for sample_id, files in samples.items():
                flattened_data.append((sample_id, f"{patient}_{sample_id}", files))

        # Sort by Snn keys
        flattened_data.sort(key=lambda x: int(x[0][1:]))

        # Construct the new dictionary
        flattened_patient_data = {f"{item[1]}": {item[0]: item[2]} for item in flattened_data}

        for patient_ref_mid, samples in flattened_patient_data.items():
            # Split the patient_ref_mid into patient_ref and MID, patient_ref can contain _ so get the last one
            patient_ref = patient_ref_mid.rsplit("_", 1)[0]
            for sample_id, data in samples.items():
                tag = data[0]
                personal_information_id, medical_information_id = patient_ids[patient_ref]
                # if a Serial number is provided, use it, otherwise do not include it in the JSON
                _bds_number = None
                if patient_ref in _bds_mapping:
                    print(f"Using Serial number {_bds_mapping[patient_ref]} for patient {patient_ref}")
                    _bds_number = _bds_mapping[patient_ref]
                    analyses.append(
                        {
                            "definition": {
                                "sampleId": sample_id,
                                "multiplexId": sample_id,
                                "sgaPipelineId": pipeline_id,
                                "userRef": patient_ref if len(tag) == 0 else f"{patient_ref}-{tag}",
                                "sampleTypeId": JsonBuilder.tag_to_sample_type_id_map_msk_access[
                                    tag] if _is_msk_pipeline else sample_type_id,
                                "libraryType": JsonBuilder.lib_type[tag],
                                "bdsNumber": _bds_number
                            },
                            "patient": {
                                "personalInformationId": personal_information_id,
                                "medicalInformationId": medical_information_id
                            },
                            "isControlSample": True if tag in ["CP", "CN"] else False,
                            "files": [
                                {
                                    "definition": {
                                        "name": file
                                    },
                                } for file in [file.replace("\\", "\\\\") for file in data[1:]]
                            ]
                        }
                    )
                else:
                    print(f"No Serial number provided for patient {patient_ref}")
                    analyses.append(
                        {
                            "definition": {
                                "sampleId": sample_id,
                                "multiplexId": sample_id,
                                "sgaPipelineId": pipeline_id,
                                "userRef": patient_ref if len(tag) == 0 else f"{patient_ref}-{tag}",
                                "sampleTypeId": JsonBuilder.tag_to_sample_type_id_map_msk_access[
                                    tag] if _is_msk_pipeline else sample_type_id,
                                "libraryType": JsonBuilder.lib_type[tag],
                            },
                            "patient": {
                                "personalInformationId": personal_information_id,
                                "medicalInformationId": medical_information_id
                            },
                            "isControlSample": True if tag in ["CP", "CN"] else False,
                            "files": [
                                {
                                    "definition": {
                                        "name": file
                                    },
                                } for file in [file.replace("\\", "\\\\") for file in data[1:]]
                            ]
                        }
                    )

        return analyses

    @staticmethod
    def _build_topologies(patient_data):
        """
        :param patient_data:
        :return: list of topologies for ADE
        """
        global debug
        topologies = []
        if debug:
            print("[DEBUG] Building topologies")
            print(json.dumps(patient_data, indent=3))
        for patient_ref, samples in patient_data.items():
            tags = [val[0] for val in samples.values()]
            tags.sort()
            # Only create a topology if there's a pair of D/R or N/T samples
            if len(tags) != 2 or (tags != ["D", "R"] and (tags != ["N", "T"] and tags != ["lib1", "lib2"])):
                continue
            topologies.append(
                {
                    "definition": {
                        "type": JsonBuilder.top_type[tags[0]]
                    },
                    "references": [
                        {
                            "analysisReference": {
                                "sampleId": sample_id
                            },
                            "role": JsonBuilder.role[data[0]],
                            "metadata": {}
                        }
                        for sample_id, data in samples.items()
                    ]
                }
            )

        return topologies


class Logger:

    @staticmethod
    def log(result):
        """
        Print the stdout and stderr of result if not empty
        :param result: return value of subprocess.run()
        :return:
        """
        out = result.stdout.strip()
        if len(out) > 0:
            print(out)

        err = result.stderr.strip()
        if len(err) > 0:
            print(err)


def confirm():
    """
    Dispplay warnings/info and return whether user wishes to continue
    :return:
    """
    print("Prerequisites")
    print(" - be logged in with sg-upload-v2-latest.jar")
    print(" - all files in target folder should use the same pipeline")
    print("Please do not upload any files containing nominative information or any other direct identifier ", end="")
    print("related to a patient (e.g. patient’s first and/or last names in file name)")

    return input("Do you wish to continue? (y/n) ").strip()[0].lower() == "y"


def run(command):
    global debug
    if debug:
        print("[DEBUG] Command to be run:", command)
    return subprocess.run(command, capture_output=True, text=True)


class RegexHelper:
    ILLUMINA_RE = r"^([^_]+?)(-(T|N|CP|CN|D|R|lib1|lib2))?_(S[0-9][0-9]*)_L\d+_R\d+_\d+\.fastq\.gz$"
    ILLUMINA_BREAKDOWN = ["0", "3", "2"]

    def __init__(self, filename):
        if filename is None:
            # Default to Illumina
            self.exprs = [(RegexHelper.ILLUMINA_RE, RegexHelper.ILLUMINA_BREAKDOWN)]
        else:
            self.exprs = RegexHelper.load_expressions(filename)

        self.functions = RegexHelper.build_functions(self.exprs)

    def update(self, filename, patients):
        updated = False
        for _func in self.functions:
            updated = _func(filename, patients)
            if updated:
                break
        return updated

    @staticmethod
    def load_expressions(file):
        exprs = []
        with open(file, "r") as file_in:
            lines = file_in.readlines()

        for _i in range(0, len(lines), 2):
            regex = lines[_i].strip()
            breakdown = lines[_i + 1].strip().split(" ")
            exprs.append((regex, breakdown))

        return exprs

    @staticmethod
    def build_functions(exprs):
        funcs = []

        for _expr in exprs:
            _regex = re.compile(_expr[0])
            _p_ref_ex = _expr[1][0]
            _mid_ex = _expr[1][1]
            _tag_ex = _expr[1][2] if len(_expr[1]) == 3 else None
            funcs.append(RegexHelper._build_function(_regex, _p_ref_ex, _mid_ex, _tag_ex))

        return funcs

    @staticmethod
    def _build_function(regex, p_ref_ex, mid_ex, tag_ex):
        def _func(filename, patients):
            basename = os.path.basename(filename)
            match = re.search(regex, basename)
            if match:
                p_ref = "".join([match.groups()[int(_c)] if _c.isdigit() else _c for _c in p_ref_ex])
                mid = "".join([match.groups()[int(_c)] if _c.isdigit() else _c for _c in mid_ex])
                tag = ""
                try:
                    tag = "".join([match.groups()[int(_c)] if _c.isdigit() else _c for _c in tag_ex])
                except TypeError:
                    # List comprehension evaluates to [None] for Illumina without a tag, so ignore this
                    pass

                PatientHelper.update_patient(patients, p_ref, mid, tag, filename)
                return True
            else:
                return False

        return _func


def validate_bds_number(bds_string):
    # strip trailing and leading whitespace
    bds_string = bds_string.strip()

    # Regex to match the required format
    pattern = r'^BDS-([0-9]{10})-([0-9]{2})$'

    match = re.match(pattern, bds_string)

    if match:
        digits = match.group(1)  # The 10-digit sequence
        last_two_digits = match.group(2)  # The last two digits
        # parse number from last two digits
        # cast to int last two digits
        last_two_digits_number = int(last_two_digits)

        # Check if the sum of the first 10 digits equals the sum of the last two digits
        sum_first_ten_digits = sum(int(digit) for digit in digits)
        if debug:
            # print all groups
            print(f"Matched groups: {match.groups()}")
            print(f"First 10 digits: {digits}")
            print(f"Last two digits: {last_two_digits}")
            print(f"Sum of the first 10 digits: {sum_first_ten_digits}")
            print(f"Last two digits number: {last_two_digits_number}")
        if sum_first_ten_digits == last_two_digits_number:
            # Sum condition met
            return True
        else:
            # Sum condition not met
            return False
    else:
        # Regex pattern not matched
        return False


if __name__ == "__main__":
    print(f"adegen.py v{VERSION}")
    # Parse the command line args
    _parser = argparse.ArgumentParser(description="Generate ADE file from FastQ folder")
    _parser.add_argument("folder", help="Path to a folder containing FastQ files")
    _parser.add_argument("-j", "--jar", default="./sg-upload-v2-latest.jar",
                         help="Location of sg-upload-v2-latest.jar (defaults to ./sg-upload-v2-latest.jar)")
    _parser.add_argument("-o", "--output", help="Output Json file (overwites without warning)")
    _parser.add_argument("-r", "--ref", help="A name for the run")
    _parser.add_argument("-p", "--pipeline", default=-1, help="ID of pipeline")
    _parser.add_argument("-s", "--sampletype", default=108000,
                         help="sampleTypeId to apply to all samples - defaults to 108000 (Peripheral Blood)")
    _parser.add_argument("-c", "--confirm", action="store_true", help="Confirm use of script")
    _parser.add_argument("-i", "--clientId", help="Client ID for the data")
    _parser.add_argument("-d", "--deep", action="store_true", help="Recurse through target folder")
    _parser.add_argument("-v", "--verbose", action="store_true", help="Debug mode")
    _parser.add_argument("-x", "--regex", help="Override regex")
    _parser.add_argument("-y", "--yaml", help="An override file for the CLI")
    _parser.add_argument("-fp", "--forceplatform", action="store_true", help="Force platform services")

    # Extended argument parsing to include Serial number mapping file and Serial number options
    _parser.add_argument("--bdsMappingFile",
                         help="Path to the mapping file containing mapping of patient references to Serial Numbers")
    _parser.add_argument("--bdsNumber", help="Mandatory Serial Number for all SOPHiA GENETICS bundle solutions")

    _parser.set_defaults(verbose=False)
    _args = _parser.parse_args()

    if _args.verbose:
        debug = True
        print("[DEBUG] Debug mode active")

    # Get user confirmation to continue
    if not _args.confirm and not confirm():
        exit(0)

    # Make sure we have the upload CLI and set the upload command
    JAR_COMMAND = ["java", "-jar"]
    if debug:  # if debug is True, add debug flag to the command
        JAR_COMMAND += ["-Dlog.level=DEBUG"]

    if os.path.exists(_args.jar):
        if _args.yaml:
            if os.path.exists(_args.yaml):
                JAR_COMMAND += [f"-Dmicronaut.config.files={_args.yaml}"]
            else:
                print(f"Couldn't find {_args.yaml}")
                exit(7)
        JAR_COMMAND += [_args.jar]
        print(JAR_COMMAND)
    else:
        print(f"Error: {_args.jar} not found")
        exit(6)

    # Get pipeline info
    print("Fetching pipeline info")
    _pipeline_id, _sequencer_id, _analysis_type_id = UserHelper.get_pipeline(JAR_COMMAND[:], _args.pipeline,
                                                                             _args.forceplatform)
    if _sequencer_id == -1:
        exit(5)

    print(f"Found Pipeline ID: {_pipeline_id}, Sequencer ID: {_sequencer_id} Analysis Type: {_analysis_type_id}")

    # Parse the patient data from files in fastqfolder
    _patient_data = PatientHelper.read(os.path.join(os.getcwd(), _args.folder), _args.deep, _args.regex)

    _is_msk_pipeline = is_msk_pipeline(_analysis_type_id)
    if _patient_data is None or len(_patient_data) == 0:
        print("No patient data found")
        exit(1)
    elif _is_msk_pipeline:
        print("MSK pipeline detected")
        if not PatientHelper.is_valid_msk(_patient_data):
            print("Patient data is not valid for MSK pipeline")
            exit(1)
    else:
        if not PatientHelper.is_valid(_patient_data):
            print("Patient data is not valid")
            exit(1)

    # Init the Serial number mapping
    bds_mapping = {}

    # If a single Serial number is provided, apply it to all samples
    if _args.bdsNumber:
        if not validate_bds_number(_args.bdsNumber):
            print("Serial number", _args.bdsNumber, "is invalid")
            sys.exit(1)
        for patient_ref, data in _patient_data.items():
            bds_mapping[patient_ref] = _args.bdsNumber
    elif _args.bdsMappingFile:
        print(f"Serial number mapping file provided: {_args.bdsMappingFile}")
        print("Reading Serial number mapping file ...")
        try:
            with open(_args.bdsMappingFile, 'r') as f:
                for line in f:
                    parts = line.strip().split(',')
                    if len(parts) == 2:
                        patient_ref, bds_number = parts
                        # validate the Serial number
                        if not validate_bds_number(bds_number):
                            print(f"Invalid Serial number found in mapping file: {bds_number}")
                            sys.exit(1)
                        bds_mapping[patient_ref] = bds_number
        except Exception as e:
            print(f"Error reading Serial number mapping file: {e}")
            sys.exit(1)

    # Check if there's a normal sample in the patient data and ask for user's confirmation
    if not PatientHelper.confirm_normal_sample_analysis(_patient_data):
        exit(1)

        # Get a dict of {patient_ref: (personalID, medicalId), . . . } for each patient
    _patients_ids = PatientHelper.create_patients(JAR_COMMAND[:], _patient_data, _args.clientId, _args.forceplatform)
    if _patients_ids is None:
        print("Error: couldn't create patients")
        exit(1)

    # Get user info
    _user_id, _client_id = UserHelper.get_user_info(JAR_COMMAND[:])
    if _args.clientId:
        _client_id = int(_args.clientId)
    if _user_id == -1 or _client_id == -1:
        exit(4)

    # Get userRef
    if _args.ref:
        _user_ref = _args.ref
    else:
        _user_ref = f'{pathlib.PurePath(_args.folder).name}_{datetime.datetime.now().strftime("%Y%m%d%H%M")}'

    if debug:
        print("Serial number provided: ", _args.bdsNumber)
        print("Serial number mapping file provided: ", _args.bdsMappingFile)
        print("Patient data:")
        print(json.dumps(_patient_data, indent=3))
        print("Serial number mapping:")
        print(json.dumps(bds_mapping, indent=3))

    # Put the ADE together
    _ade_json = JsonBuilder.build_json(_user_ref, _user_id, _client_id,
                                       _pipeline_id, _sequencer_id, _patient_data, _patients_ids, _args.sampletype,
                                       _is_msk_pipeline, bds_mapping)

    # Write to file if given, otherwise print to console
    if _args.output:
        save_path = os.path.join(os.getcwd(), _args.output)
        with open(save_path, "w") as f_out:
            f_out.write(_ade_json)
        print(f"Json written to {save_path}")
    else:
        print(_ade_json)
