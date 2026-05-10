.data
input_file:   .asciiz "input.txt"
desired_file: .asciiz "desired.txt"
output_file:  .asciiz "output.txt"
newline:      .asciiz "\n"

buf_size:     .word 32768
buffer:       .space 32768

NUM_SAMPLES:  .word 10

desired:      .space 40      # 10 floats
input:        .space 40      # 10 floats
crosscorr:    .space 40
autocorr:     .space 40
R:            .space 400     # 10x10 float
coeff:        .space 40
output:       .space 40

mmse:         .float 0.0

zero_f:       .float 0.0
one_f:        .float 1.0
ten:          .float 10.0
hundred:      .float 100.0
half:         .float 0.5
minus_half:   .float -0.5
zero:         .float 0.0

header_filtered: .asciiz "Filtered output: "
header_mmse:     .asciiz "\nMMSE: "
space_str:       .asciiz " "

str_buf:    .space 32
temp_str:   .space 32

error_open: .asciiz "Error: Can not open file\n"
error_read: .asciiz "Error: Can not read the file\n"
error_size: .asciiz "Error: size not match\n"

.text
.globl main

main:
    ## TODO CODE
    # --- BƯỚC 1: LOAD DỮ LIỆU ---

    la   $a0, input_file
    la   $a1, input
    jal  read_and_parse_file
    move $s7, $v0           # Lưu số lượng mẫu đọc được từ input

    # Đọc file desired.txt
    la   $a0, desired_file
    la   $a1, desired
    jal  read_and_parse_file
    
    # Kiểm tra kích thước 
    lw   $t0, NUM_SAMPLES
    bne  $s7, $t0, size_error
    bne  $v0, $t0, size_error
    j    start_processing

size_error:
    jal  write_error_to_file
    li   $v0, 17            
    li   $a0, 0            
    syscall

start_processing:

    # --- BƯỚC 2: TÍNH TOÁN TƯƠNG QUAN ---
    la   $a0, desired
    la   $a1, input
    la   $a2, crosscorr
    lw   $a3, NUM_SAMPLES
    jal  computeCrosscorrelation

    la   $a0, input
    la   $a1, autocorr
    lw   $a3, NUM_SAMPLES
    jal  computeAutocorrelation

    # --- BƯỚC 3: TẠO MA TRẬN TOEPLITZ ---
    la   $a0, autocorr
    la   $a1, R
    lw   $a2, NUM_SAMPLES
    jal  createToeplitzMatrix

    # --- BƯỚC 4: GIẢI HỆ PHƯƠNG TRÌNH & LỌC ---
    la   $a0, R
    la   $a1, crosscorr
    la   $a2, coeff
    lw   $a3, NUM_SAMPLES
    jal  solveLinearSystem

    la   $a0, input
    la   $a1, coeff
    la   $a2, output
    lw   $a3, NUM_SAMPLES
    jal  applyWienerFilter

    # --- BƯỚC 5: TÍNH MMSE ---
    jal  calculate_mmse_val

    # --- BƯỚC 6: XUẤT FILE ---
    jal  write_to_output_file

    # Thoát chương trình
    li   $v0, 10
    syscall

computeAutocorrelation:
    move $t0, $a0
    move $t1, $a1
    move $t2, $a3
    li   $t3, 0

autocorr_loop:
    beq  $t3, $t2, autocorr_done
    la   $t8, zero_f
    lwc1 $f4, 0($t8)
    sub  $t4, $t2, $t3
    li   $t5, 0

autocorr_inner:
    beq  $t5, $t4, autocorr_inner_done
    sll  $t6, $t5, 2
    addu $t7, $t0, $t6
    lwc1 $f6, 0($t7)
    addu $t9, $t5, $t3
    sll  $t9, $t9, 2
    addu $t9, $t0, $t9
    lwc1 $f8, 0($t9)
    mul.s $f10, $f6, $f8
    add.s $f4, $f4, $f10
    addi $t5, $t5, 1
    j    autocorr_inner

autocorr_inner_done:
    mtc1 $t2, $f2
    cvt.s.w $f2, $f2
    div.s $f12, $f4, $f2
    sll  $t7, $t3, 2
    addu $t7, $t1, $t7
    swc1 $f12, 0($t7)
    addi $t3, $t3, 1
    j    autocorr_loop

autocorr_done:
    jr $ra

computeCrosscorrelation:
    move $t0, $a0
    move $t1, $a1
    move $t2, $a2
    move $t3, $a3
    li   $t4, 0

cross_loop:
    beq  $t4, $t3, cross_done
    la   $t9, zero_f
    lwc1 $f4, 0($t9)
    sub  $t5, $t3, $t4
    li   $t6, 0

cross_inner:
    beq  $t6, $t5, cross_inner_done
    addu $t7, $t6, $t4
    sll  $t8, $t7, 2
    addu $t8, $t0, $t8
    lwc1 $f6, 0($t8)
    sll  $t9, $t6, 2
    addu $t9, $t1, $t9
    lwc1 $f8, 0($t9)
    mul.s $f10, $f6, $f8
    add.s $f4, $f4, $f10
    addi $t6, $t6, 1
    j    cross_inner

cross_inner_done:
    mtc1 $t3, $f2
    cvt.s.w $f2, $f2
    div.s $f12, $f4, $f2
    sll  $t8, $t4, 2
    addu $t8, $t2, $t8
    swc1 $f12, 0($t8)
    addi $t4, $t4, 1
    j    cross_loop
cross_done:
    jr $ra

createToeplitzMatrix:
    move $t0, $a0
    move $t1, $a1
    move $t2, $a2
    li   $t3, 0

toe_outer:
    beq  $t3, $t2, toe_done
    li   $t4, 0

toe_inner:
    beq  $t4, $t2, toe_next_i
    sub  $t5, $t3, $t4
    abs  $t5, $t5
    sll  $t6, $t5, 2
    addu $t6, $t0, $t6
    lwc1 $f6, 0($t6)
    mul  $t7, $t3, $t2
    add  $t7, $t7, $t4
    sll  $t7, $t7, 2
    addu $t7, $t1, $t7
    swc1 $f6, 0($t7)
    addi $t4, $t4, 1
    j    toe_inner

toe_next_i:
    addi $t3, $t3, 1
    j    toe_outer

toe_done:
    jr $ra

solveLinearSystem:
    # Backup s-registers
    addi $sp, $sp, -16
    sw   $s0, 0($sp)
    sw   $s1, 4($sp)
    sw   $s2, 8($sp)
    sw   $s3, 12($sp)
    
    move $s0, $a0
    move $s1, $a1
    move $s2, $a2
    move $s3, $a3
    li   $t0, 0

elim_outer:
    beq  $t0, $s3, elim_done
    mul  $t2, $t0, $s3
    add  $t2, $t2, $t0
    sll  $t2, $t2, 2
    addu $t2, $s0, $t2
    lwc1 $f16, 0($t2)        # pivot
    addi $t3, $t0, 1

elim_k:
    beq  $t3, $s3, elim_next_i
    mul  $t4, $t3, $s3
    add  $t4, $t4, $t0
    sll  $t4, $t4, 2
    addu $t4, $s0, $t4
    lwc1 $f18, 0($t4)
    div.s $f20, $f18, $f16   # factor
    move $t5, $t0

elim_j:
    beq  $t5, $s3, elim_b
    mul  $t6, $t3, $s3
    add  $t6, $t6, $t5
    sll  $t6, $t6, 2
    addu $t6, $s0, $t6
    lwc1 $f22, 0($t6)
    mul  $t7, $t0, $s3
    add  $t7, $t7, $t5
    sll  $t7, $t7, 2
    addu $t7, $s0, $t7
    lwc1 $f24, 0($t7)
    mul.s $f26, $f20, $f24
    sub.s $f22, $f22, $f26
    swc1 $f22, 0($t6)
    addi $t5, $t5, 1
    j    elim_j

elim_b:
    sll  $t8, $t3, 2
    addu $t8, $s1, $t8
    lwc1 $f28, 0($t8)
    sll  $t9, $t0, 2
    addu $t9, $s1, $t9
    lwc1 $f30, 0($t9)
    mul.s $f26, $f20, $f30
    sub.s $f28, $f28, $f26
    swc1 $f28, 0($t8)
    addi $t3, $t3, 1
    j    elim_k

elim_next_i:
    addi $t0, $t0, 1
    j    elim_outer

elim_done:
    # Back substitution
    addi $t0, $s3, -1

back_i:
    bltz $t0, solve_final
    sll  $t1, $t0, 2
    addu $t1, $s1, $t1
    lwc1 $f20, 0($t1)
    addi $t2, $t0, 1

back_j:
    bge  $t2, $s3, back_div
    mul  $t3, $t0, $s3
    add  $t3, $t3, $t2
    sll  $t3, $t3, 2
    addu $t3, $s0, $t3
    lwc1 $f22, 0($t3)
    sll  $t4, $t2, 2
    addu $t4, $s2, $t4
    lwc1 $f24, 0($t4)
    mul.s $f26, $f22, $f24
    sub.s $f20, $f20, $f26
    addi $t2, $t2, 1
    j    back_j

back_div:
    mul  $t5, $t0, $s3
    add  $t5, $t5, $t0
    sll  $t5, $t5, 2
    addu $t5, $s0, $t5
    lwc1 $f22, 0($t5)
    div.s $f20, $f20, $f22
    sll  $t6, $t0, 2
    addu $t6, $s2, $t6
    swc1 $f20, 0($t6)
    addi $t0, $t0, -1
    j    back_i

solve_final:
    lw   $s0, 0($sp)
    lw   $s1, 4($sp)
    lw   $s2, 8($sp)
    lw   $s3, 12($sp)
    addi $sp, $sp, 16
    jr   $ra

applyWienerFilter:
    move $t0, $a0
    move $t1, $a1
    move $t2, $a2
    move $t3, $a3
    li   $t4, 0

apply_outer:
    beq  $t4, $t3, apply_done
    la   $t9, zero_f
    lwc1 $f20, 0($t9)
    li   $t5, 0

apply_inner:
    bgt  $t5, $t4, apply_store
    sll  $t6, $t5, 2
    addu $t6, $t1, $t6
    lwc1 $f22, 0($t6)
    sub  $t7, $t4, $t5
    sll  $t7, $t7, 2
    addu $t7, $t0, $t7
    lwc1 $f24, 0($t7)
    mul.s $f26, $f22, $f24
    add.s $f20, $f20, $f26
    addi $t5, $t5, 1
    j    apply_inner

apply_store:
    sll  $t8, $t4, 2
    addu $t8, $t2, $t8
    swc1 $f20, 0($t8)
    addi $t4, $t4, 1
    j    apply_outer

apply_done:
    jr $ra

# =========================================================
# CÁC HÀM (Xuất dữ liệu và Chuyển đổi)
# =========================================================

write_to_output_file:
    # Backup ra
    addi $sp, $sp, -4
    sw   $ra, 0($sp)

    li   $v0, 13
    la   $a0, output_file
    li   $a1, 1
    li   $a2, 0
    syscall
    move $s0, $v0        # File descriptor

    # Write Header Filtered
    li   $v0, 15
    move $a0, $s0
    la   $a1, header_filtered
    li   $a2, 17
    syscall

    # Loop write output
    lw   $s1, NUM_SAMPLES
    li   $s2, 0

write_loop:
    beq  $s2, $s1, write_mmse_part
    
    # 1. Lấy dữ liệu và chuyển đổi float thành chuỗi
    sll  $t0, $s2, 2
    la   $t1, output
    addu $t1, $t1, $t0
    lwc1 $f12, 0($t1)
    la   $a0, str_buf
    li   $a1, 1              # Giữ 1 chữ số thập phân 
    jal  float_to_str
    
    # 2. Ghi số đã chuyển đổi vào file
    move $t8, $v0            # Lưu độ dài chuỗi số
    li   $v0, 15
    move $a0, $s0
    la   $a1, str_buf
    move $a2, $t8
    syscall

    # 3. Kiểm tra để không ghi dấu cách ở phần tử cuối
    addi $t9, $s1, -1        # t9 = NUM_SAMPLES - 1
    beq  $s2, $t9, skip_space # Nếu đang ở phần tử cuối cùng, nhảy qua đoạn ghi dấu cách
    
    li   $v0, 15
    move $a0, $s0
    la   $a1, space_str
    li   $a2, 1
    syscall

skip_space:
    addi $s2, $s2, 1
    j    write_loop

write_mmse_part:
    # 4. Ghi tiêu đề MMSE 
    li   $v0, 15
    move $a0, $s0
    la   $a1, header_mmse
    li   $a2, 7              # Chuỗi "\nMMSE: " có đúng 7 ký tự
    syscall

    # 5. Ghi giá trị MMSE (1 chữ số thâp phân)
    lwc1 $f12, mmse
    la   $a0, str_buf
    li   $a1, 1             
    jal  float_to_str
    
    move $t8, $v0
    li   $v0, 15
    move $a0, $s0
    la   $a1, str_buf
    move $a2, $t8
    syscall

    # 6. Đóng file và kết thúc
    li   $v0, 16
    move $a0, $s0
    syscall

    lw   $ra, 0($sp)
    addi $sp, $sp, 4
    jr   $ra

# ---------------------------------------------------------
# float_to_str(buffer $a0, value $f12, decimals $a1) -> length $v0
# ---------------------------------------------------------
float_to_str:
    # --- PROLOGUE: Backup registers to Stack ---
    addi $sp, $sp, -40      # Tạo khoảng trống cho 10 thanh ghi (40 bytes)
    sw   $ra, 36($sp)
    sw   $s0, 32($sp)
    sw   $s1, 28($sp)
    sw   $s2, 24($sp)
    sw   $s3, 20($sp)
    sw   $s4, 16($sp)
    sw   $s5, 12($sp)
    sw   $s6, 8($sp)
    sw   $s7, 4($sp)
    sw   $t0, 0($sp)

    # Chuyển dữ liệu vào các thanh ghi an toàn
    move $s0, $a0           # $s0 = buffer ptr
    move $s6, $a1           # $s6 = decimals (1 or 2)

    # Chọn multiplier
    li   $t0, 1
    beq  $s6, $t0, ft_one
    # decimals == 2
    la   $t1, hundred
    lwc1 $f1, 0($t1)
    j    ft_after_mul
ft_one:
    la   $t1, ten
    lwc1 $f1, 0($t1)
ft_after_mul:
    mul.s $f2, $f12, $f1    # f2 = value * mul

    # Làm tròn (Rounding)
    la   $t2, zero_f
    lwc1 $f0, 0($t2)
    c.lt.s $f2, $f0
    bc1t ft_neg
    la   $t3, half
    lwc1 $f3, 0($t3)
    add.s $f2, $f2, $f3
    j    ft_after_round
ft_neg:
    la   $t3, half
    lwc1 $f3, 0($t3)
    sub.s $f2, $f2, $f3

ft_after_round:
    cvt.w.s $f4, $f2
    mfc1 $s1, $f4           # $s1 = integer scaled value (thay $t4)

    # Xử lý dấu
    li   $s2, 0             # $s2 = negative flag (thay $t5)
    bltz $s1, ft_neg_int
    j    ft_pos_int
ft_neg_int:
    li   $s2, 1
    subu $s1, $zero, $s1
ft_pos_int:
    # Xác định số chia (tương ứng decimals)
    li   $t6, 10
    li   $t0, 1
    beq  $s6, $t0, ft_div_ready
    li   $t6, 100
ft_div_ready:
    div  $s1, $t6
    mflo $s3                # $s3 = integer_part (thay $t8)
    mfhi $s4                # $s4 = remainder/frac_part (thay $t9)

    move $s5, $s0           # $s5 = current buffer ptr (thay $t10)
    li   $v0, 0             # length counter

    # Ghi dấu '-' nếu là số âm
    beq  $s2, $zero, ft_write_int
    li   $t7, 45            # ASCII '-'
    sb   $t7, 0($s5)
    addi $s5, $s5, 1
    addi $v0, $v0, 1

ft_write_int:
    # convert integer_part to string
    move $t1, $s3           # integer_part
    li   $t2, 0             # digit counter
    la   $t3, temp_str
    bnez $t1, ft_int_digits

    # Trường hợp integer_part = 0
    li   $t7, 48            # ASCII '0'
    sb   $t7, 0($t3)
    addi $t2, $t2, 1       
    j    ft_flush_int

ft_int_digits:
    li   $t7, 10
    div  $t1, $t7
    mflo $t1
    mfhi $t4
    addi $t4, $t4, 48       # Chuyển sang mã ASCII
    sb   $t4, 0($t3)
    addi $t3, $t3, 1
    addi $t2, $t2, 1        # Tăng biến đếm chữ số
    bnez $t1, ft_int_digits

ft_flush_int:
    la   $t3, temp_str
    add  $t3, $t3, $t2
    addi $t3, $t3, -1

ft_rev_loop:
    blez $t2, ft_after_int
    lb   $t4, 0($t3)
    sb   $t4, 0($s5)
    addi $s5, $s5, 1
    addi $v0, $v0, 1
    addi $t3, $t3, -1
    addi $t2, $t2, -1
    j    ft_rev_loop

ft_after_int:
    li   $t7, 46            # '.'
    sb   $t7, 0($s5)
    addi $s5, $s5, 1
    addi $v0, $v0, 1

    # Xử lý phần thập phân
    li   $t0, 1
    beq  $s6, $t0, ft_frac_one_digit
    # Decimals == 2
    li   $t7, 10
    div  $s4, $t7
    mflo $t1                # Digit 1
    mfhi $t2                # Digit 2
    addi $t1, $t1, 48
    sb   $t1, 0($s5)
    addi $t2, $t2, 48
    sb   $t2, 1($s5)
    addi $s5, $s5, 2
    addi $v0, $v0, 2
    j    ft_done

ft_frac_one_digit:
    addi $t1, $s4, 48
    sb   $t1, 0($s5)
    addi $s5, $s5, 1
    addi $v0, $v0, 1

ft_done:
    # --- EPILOGUE: Restore registers ---
    lw   $t0, 0($sp)
    lw   $s7, 4($sp)
    lw   $s6, 8($sp)
    lw   $s5, 12($sp)
    lw   $s4, 16($sp)
    lw   $s3, 20($sp)
    lw   $s2, 24($sp)
    lw   $s1, 28($sp)
    lw   $s0, 32($sp)
    lw   $ra, 36($sp)
    addi $sp, $sp, 40
    jr   $ra

calculate_mmse_val:
    la   $t0, desired
    la   $t1, output
    lw   $t2, NUM_SAMPLES
    la   $t3, zero_f
    lwc1 $f0, 0($t3)
    li   $t4, 0
mmse_l:
    beq  $t4, $t2, mmse_d
    sll  $t5, $t4, 2
    addu $t6, $t0, $t5
    lwc1 $f1, 0($t6)
    addu $t7, $t1, $t5
    lwc1 $f2, 0($t7)
    sub.s $f3, $f1, $f2
    mul.s $f3, $f3, $f3
    add.s $f0, $f0, $f3
    addi $t4, $t4, 1
    j    mmse_l
mmse_d:
    mtc1 $t2, $f4
    cvt.s.w $f4, $f4
    div.s $f0, $f0, $f4
    swc1 $f0, mmse
    jr   $ra

# ---------------------------------------------------------
# read_and_parse_file(filename $a0, array_ptr $a1) -> count $v0
# ---------------------------------------------------------
read_and_parse_file:
    addi $sp, $sp, -20
    sw   $ra, 16($sp)
    sw   $s0, 12($sp)
    sw   $s1, 8($sp)
    sw   $s2, 4($sp)
    
    move $s0, $a1           # Array pointer
    # Open file
    li   $v0, 13
    li   $a1, 0             # Read-only
    li   $a2, 0
    syscall
    move $s1, $v0           # File descriptor
    
    # Read entire file to buffer
    li   $v0, 14
    move $a0, $s1
    la   $a1, buffer
    lw   $a2, buf_size
    syscall
    move $s2, $v0           # Number of bytes read
    
    # Close file
    li   $v0, 16
    move $a0, $s1
    syscall
    
    # Parse buffer to floats
    la   $a0, buffer
    move $a1, $s0
    move $a2, $s2           # Buffer length
    jal  parse_ascii_to_float
    
    lw   $s2, 4($sp)
    lw   $s1, 8($sp)
    lw   $s0, 12($sp)
    lw   $ra, 16($sp)
    addi $sp, $sp, 20
    jr   $ra

# Hàm xử lý chuỗi ký tự thành số thực 
parse_ascii_to_float:
    # parse_ascii_to_float(buffer_ptr $a0, array_ptr $a1, len $a2)
    # Returns count in $v0
    addi $sp, $sp, -32
    sw   $ra, 28($sp)
    sw   $s0, 24($sp)
    sw   $s1, 20($sp)
    sw   $s2, 16($sp)
    sw   $s3, 12($sp)
    sw   $s4, 8($sp)
    sw   $s5, 4($sp)

    move $s0, $a0      # buffer ptr
    move $s1, $a1      # array ptr
    move $s2, $a2      # buffer len
    li   $v0, 0        # count

    li   $s3, 0        # index into buffer
parse_loop:
    bge  $s3, $s2, parse_done
    addu $t0, $s0, $s3
    lb   $t1, 0($t0)
    # skip whitespace
    li   $t2, 32
    beq  $t1, $t2, parse_skip_ws
    li   $t2, 9
    beq  $t1, $t2, parse_skip_ws
    li   $t2, 10
    beq  $t1, $t2, parse_skip_ws
    li   $t2, 13
    beq  $t1, $t2, parse_skip_ws
    # start parsing a number token
    li   $t3, 0        # sign flag 0=+,1=-
    li   $t4, 0        # integer_accum
    li   $t5, 0        # frac_accum
    li   $t6, 0        # frac_count
    # handle sign
    li   $t2, 45       # '-'
    beq  $t1, $t2, parse_sign
    j    parse_digits_loop
parse_sign:
    li   $t3, 1
    addi $s3, $s3, 1
    beq  $s3, $s2, parse_done
    addu $t0, $s0, $s3
    lb   $t1, 0($t0)

parse_digits_loop:
    bge  $s3, $s2, parse_finish_token
    # if digit
    li   $t2, 48
    li   $t7, 57
    blt  $t1, $t2, parse_after_digit_check
    bgt  $t1, $t7, parse_after_digit_check
    li   $t9, 10
    mul  $t4, $t4, $t9
    addi $t9, $t1, -48
    add  $t4, $t4, $t9
    addi $s3, $s3, 1
    beq  $s3, $s2, parse_finish_token
    addu $t0, $s0, $s3
    lb   $t1, 0($t0)
    j    parse_digits_loop

parse_after_digit_check:
    # check decimal point
    li   $t2, 46
    beq  $t1, $t2, parse_frac_start
    # token end
    j    parse_finish_token

parse_frac_start:
    addi $s3, $s3, 1
    beq  $s3, $s2, parse_finish_token
    addu $t0, $s0, $s3
    lb   $t1, 0($t0)
parse_frac_loop:
    li   $t2, 48
    li   $t7, 57
    blt  $t1, $t2, parse_finish_token
    bgt  $t1, $t7, parse_finish_token
    li   $t9, 10
    mul  $t5, $t5, $t9
    addi $t9, $t1, -48
    add  $t5, $t5, $t9
    addi $t6, $t6, 1
    addi $s3, $s3, 1
    beq  $s3, $s2, parse_finish_token
    addu $t0, $s0, $s3
    lb   $t1, 0($t0)
    j    parse_frac_loop

parse_finish_token:
    # convert integer_accum and frac_accum to float in $f0
    # f = integer_accum + frac_accum / 10^frac_count
    mtc1 $t4, $f2
    cvt.s.w $f2, $f2
    beq  $t6, $zero, parse_store_float
    mtc1 $t5, $f4
    cvt.s.w $f4, $f4
    # compute pow10
    li   $t7, 1
    li   $t8, 0
pow10_loop:
    beq  $t8, $t6, pow10_done
    li   $t9, 10
    mul  $t7, $t7, $t9
    addi $t8, $t8, 1
    j    pow10_loop
pow10_done:
    mtc1 $t7, $f6
    cvt.s.w $f6, $f6
    div.s $f4, $f4, $f6
    add.s $f2, $f2, $f4

parse_store_float:
    # apply sign
    beq  $t3, $zero, parse_skip_sign
    neg.s $f2, $f2
parse_skip_sign:
    # store float into array
    sll  $t0, $v0, 2
    addu $t0, $s1, $t0
    swc1 $f2, 0($t0)
    addi $v0, $v0, 1
    # advance past token delimiter if not at end
    j    parse_loop

parse_skip_ws:
    addi $s3, $s3, 1
    j    parse_loop

parse_done:
    lw   $ra, 28($sp)
    lw   $s0, 24($sp)
    lw   $s1, 20($sp)
    lw   $s2, 16($sp)
    lw   $s3, 12($sp)
    lw   $s4, 8($sp)
    lw   $s5, 4($sp)
    addi $sp, $sp, 32
    jr   $ra

write_error_to_file:
    li   $v0, 13
    la   $a0, output_file
    li   $a1, 1
    syscall
    move $t0, $v0
    
    li   $v0, 15
    move $a0, $t0
    la   $a1, error_size
    li   $a2, 22            # 22 (đủ cho chuỗi "Error: size not match\n")
    syscall
    
    li   $v0, 16
    move $a0, $t0
    syscall
    jr   $ra

