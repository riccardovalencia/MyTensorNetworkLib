# legacy

Old code kept for reference only. It is **not** compiled into the library and some files do not compile
(they depend on an external TDVP implementation or on removed declarations).

- `spin_boson/TEBD_backup.*`: old monolithic copy of the TEBD routines, now split in `spin_boson/TEBD*`.
- `spin_boson/MyTDVP_testing.*`: TDVP-based Lindblad evolution (needs the ITensor TDVP add-on).
- `spin_boson/kondo_model.*`: unfinished Kondo-model code (duplicates the classes in `core/MyClasses`).
- `spin_boson/TEBD_testing.*`: test variants (`*_testing`) of the Rydberg and Lindblad gate builders.
