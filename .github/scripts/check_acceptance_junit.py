"""Fail CI if a critical module is absent or its success-path tests skip."""
import sys
import xml.etree.ElementTree as ET

REQUIRED = {
    "tests.test_import_circ", "tests.test_import_ciri3_native", "tests.test_qc_report",
    "tests.test_mudata_export", "tests.test_mudata_downstream",
    "tests.test_mudata_sync_semantics", "tests.test_export_multimodal",
    "tests.test_full_length_workflow", "tests.test_smartseq3_ciri3_workflow",
}
# This negative dependency test applies only when MuData is absent.
ALLOWED_SKIPS = {
    ("tests.test_mudata_export", "test_export_mudata_missing_package_gives_clear_error"),
    ("tests.test_smartseq3_ciri3_workflow", "test_mudata_requested_without_package_raises"),
}


def main(path):
    cases = list(ET.parse(path).iter("testcase"))
    for module in sorted(REQUIRED):
        selected = [c for c in cases if c.get("classname") == module]
        if not selected:
            raise SystemExit(f"Missing acceptance module: {module}")
        for case in selected:
            skipped = case.find("skipped")
            identity = (module, case.get("name"))
            if skipped is not None and not (identity in ALLOWED_SKIPS and skipped.get("message") == "mudata installed"):
                raise SystemExit(f"Unexpected acceptance skip: {identity}: {skipped.attrib}")
            if case.find("failure") is not None or case.find("error") is not None:
                raise SystemExit(f"Acceptance failure: {identity}")
        skipped_count = sum(c.find("skipped") is not None for c in selected)
        print(f"{module}: {len(selected) - skipped_count} executed, {skipped_count} explicitly allowed skips")


if __name__ == "__main__":
    main(sys.argv[1])
