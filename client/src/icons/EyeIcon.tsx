import Icon from './svg/eye.svg';

function EyeIcon({ size }: { size: number }) {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="eye icon" />
    </div>
  );
}

export default EyeIcon;
