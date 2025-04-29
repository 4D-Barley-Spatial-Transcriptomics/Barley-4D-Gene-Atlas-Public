session_info <- sessionInfo()
session_info_str <- capture.output(print(session_info))

file_path <- "FELIPE_session_info.txt"
writeLines(session_info_str, file_path)

cat("Session information has been written to", file_path)
