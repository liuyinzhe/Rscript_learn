#!/usr/bin/env Rscript

# 获取所有命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数是否存在
if (length(args) == 0) {
  stop("错误：未提供路径参数！用法：Rscript script.R <路径>")
}

# 提取第一个参数作为路径
input_path <- args[1]

# 验证路径是否存在
if (!file.exists(input_path)) {
  stop(paste("错误：路径", input_path, "不存在！"))
}

# 输出结果
cat("成功读取路径：", input_path, "\n")
cat("路径类型：", ifelse(dir.exists(input_path), "目录", "文件"), "\n")


##############################################################################

#!/usr/bin/env Rscript

# ---- 获取脚本路径 ----
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
  if (length(script_path) == 0) stop("请使用Rscript运行此脚本")
  return(normalizePath(script_path))
}

# ---- 主逻辑 ----
main <- function() {
  # 获取脚本目录
  script_dir <- dirname(get_script_path())
  
  # 拼接数据文件路径
  data_file <- file.path(script_dir, "input_data.csv")
  cat("正在读取数据文件：", data_file, "\n")
  
  # 验证文件存在性
  if (!file.exists(data_file)) {
    stop(paste("数据文件不存在：", data_file))
  }
  
  # 读取数据
  data <- read.csv(data_file)
  cat("成功读取", nrow(data), "行数据\n")
  
  # 处理数据（示例：计算均值）
  avg_value <- mean(data$value, na.rm = TRUE)
  cat("平均值：", avg_value, "\n")
  
  # 输出结果到同目录
  output_file <- file.path(script_dir, "result.txt")
  writeLines(paste("平均值：", avg_value), output_file)
  cat("结果已保存至：", output_file, "\n")
}

# 执行主函数
main()
