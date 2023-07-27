resource "aws_iam_role" "batch_execution_role" {
  name = "${data.aws_caller_identity.current.account_id}-${var.batch_import_role_name}"

  assume_role_policy = jsonencode({
    Version   = "2012-10-17"
    Statement = [
      {
        Effect    = "Allow",
        Principal = {
          Service = "batch.amazonaws.com"
        },
        Action = "sts:AssumeRole"
      },
      {
        Effect    = "Allow",
        Principal = {
          Service = "ecs-tasks.amazonaws.com"
        },
        Action = "sts:AssumeRole"
      }
    ]
  })

  inline_policy {
    name = "${data.aws_caller_identity.current.account_id}-${var.batch_import_policy_name}"

    policy = jsonencode({
      Version   = "2012-10-17"
      Statement = [
        {
          "Effect" = "Allow",
          "Action" = [
            "s3:GetObject",
            "s3:ListBucket"
          ],
          "Resource" = [
            var.bucket_arn,
          ]
        },
        {
          "Effect" = "Allow",
          "Action" = [
            "secretsmanager:GetSecretValue",
            "secretsmanager:DescribeSecret",
            "secretsmanager:ListSecretVersionIds"
          ],
          "Resource" = [var.db_secret_arn]
        },
        {
          "Effect" = "Allow",
          "Action" = [
            "logs:DescribeLogGroups"
          ],
          "Resource" = ["*"]
        }
      ]
    })
  }
}

resource "aws_iam_role_policy_attachment" "batch_role_policy_attachment_ecr" {
  role       = aws_iam_role.batch_execution_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryReadOnly"
}

resource "aws_iam_role_policy_attachment" "batch_role_policy_attachment_batch_sr" {
  role       = aws_iam_role.batch_execution_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
}


resource "aws_iam_role_policy_attachment" "batch_role_policy_attachment_ecs_task" {
  role       = aws_iam_role.batch_execution_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}


resource "aws_iam_role_policy" "manage_ecs" {
  name = "${data.aws_caller_identity.current.account_id}-${var.manage_ecs_policy_name}"
  role = aws_iam_role.batch_execution_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = [
          "ecs:DeleteCluster", "ecs:DescribeClusters", "ecs:ListClusters", "ecs:DeleteCluster"
        ]
        Effect   = "Allow"
        Resource = [aws_batch_compute_environment.cellenium_study_import_compute_environment.ecs_cluster_arn]
      },
    ]
  })
}