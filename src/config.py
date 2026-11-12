from pathlib import Path
import gzip
import zipfile


class FileTool:
    @staticmethod
    def get_project_dir(sub=None, *, relative=False):
        """
        Root directory of the PyGWater repository (the folder that contains src/, examples/, studies/, installation/).
        Resolved from the location of this file, so it no longer depends on the current working directory.
        """
        dir_of_project = Path(__file__).resolve().parents[1]

        if relative:
            import os
            result = Path(os.path.relpath(dir_of_project, Path().absolute()))
        else:
            result = dir_of_project

        if sub is not None:
            result /= sub

        return result

    CONFIG_FILE = 'installation/pygwater_paths.cfg'

    @staticmethod
    def _setting(key, env):
        """
        The data-root setting, in order of precedence:
          1. environment variable <env>
          2. '<key> = <path>' in <project>/installation/pygwater_paths.cfg (local, not in git; template:
             installation/pygwater_paths.cfg.example)
        :return: Path or None
        """
        import os
        value = os.environ.get(env)
        if value:
            return Path(value).expanduser()
        cfg = FileTool.get_project_dir(FileTool.CONFIG_FILE)
        if cfg.is_file():
            for line in cfg.read_text().splitlines():
                line = line.split('#', 1)[0].strip()
                if '=' in line:
                    k, v = (s.strip() for s in line.split('=', 1))
                    if k == key and v:
                        return Path(v).expanduser()
        return None

    @staticmethod
    def get_data_dir(sub=None):
        """
        External data root, kept outside the code (cf. PyGLDA_v2_external_data). ALL data PyGWater reads or writes
        live under it; nothing is written inside the code folder. Set it with 'data_root = ...' in
        installation/pygwater_paths.cfg or the environment variable PYGWATER_DATA (which takes precedence). It must
        lie outside the code folder; without a setting an error is raised (there is no default inside the code).
        Layout, by role:
            input/    external sources, only read by PyGWater: LLN/, land_mask/, mascon/, shp/, stations/, reference/
            cache/    built by PyGWater, regenerable and reused: LGD_template*/, mask/, shp/ (generated grids)
            output/   results, one folder per study (paper_JGR2025/, inversion_China/, ...)
            temp/     scratch, safe to delete
        """
        result = FileTool._setting('data_root', 'PYGWATER_DATA')
        if result is None:
            raise FileNotFoundError(
                "PyGWater data root not set. Code and data are kept apart: set 'data_root = <path>' in %s "
                "(copy installation/pygwater_paths.cfg.example) or the environment variable PYGWATER_DATA."
                % FileTool.get_project_dir(FileTool.CONFIG_FILE))
        if FileTool.get_project_dir() in [result, *result.parents]:
            raise ValueError('The data root %s lies inside the code folder %s; put it outside.'
                             % (result, FileTool.get_project_dir()))
        if sub is not None:
            result /= sub
        return result

    @staticmethod
    def _under(folder, sub):
        d = FileTool.get_data_dir(folder)
        return d if sub is None else d / sub

    @staticmethod
    def get_input_dir(sub=None):
        """external sources, only read: <data root>/input (LLN/, land_mask/, mascon/, shp/, stations/, reference/)"""
        return FileTool._under('input', sub)

    @staticmethod
    def get_mascon_dir(sub=None):
        """GRACE mascons: NetCDF (e.g. CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc, may be a symbolic link)
        and the 2-degree text extract mascon.txt (lat, lon, EWH [cm]): <data root>/input/mascon"""
        return FileTool.get_input_dir('mascon') if sub is None else FileTool.get_input_dir('mascon') / sub

    @staticmethod
    def get_cache_dir(sub=None):
        """built by PyGWater, regenerable and reused: <data root>/cache (LGD_template*/, mask/, shp/)"""
        return FileTool._under('cache', sub)

    @staticmethod
    def get_output_dir(sub=None):
        """results, one folder per study: <data root>/output"""
        return FileTool._under('output', sub)

    @staticmethod
    def get_temp_dir(sub=None):
        """scratch files, safe to delete: <data root>/temp"""
        return FileTool._under('temp', sub)

    @staticmethod
    def get_files_in_dir(fp, sub=False):
        assert fp.is_dir, f"{fp} is not a directory."

        file_list = []

        iterlist = list(fp.iterdir())
        for i in range(len(iterlist)):
            if iterlist[i].is_file():
                file_list.append(iterlist[i])

            elif sub:
                file_list += FileTool.get_files_in_dir(iterlist[i], sub=sub)

        return file_list

    @staticmethod
    def un_gz(gz_file_path: Path, target_path: Path = None):

        if target_path is None:
            target_path = gz_file_path.parent / gz_file_path.name.replace('.gz', '')

        g_file = gzip.GzipFile(gz_file_path)
        open(target_path, "wb+").write(g_file.read())
        g_file.close()

    @staticmethod
    def un_zip(zip_file_path: Path, target_path: Path = None):
        if target_path is None:
            target_path = zip_file_path.parent / zip_file_path.name.replace('.zip', '')

        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(target_path)
